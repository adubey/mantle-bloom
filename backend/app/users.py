"""Per-user store behind #560's paid tier: users, linked identities, the subscription
cache, and world ownership (#620).

Two backends share the `UserStore` interface. `SqliteUserStore` is the durable one:
a SQLite file chosen by `MANTLE_BLOOM_USER_DB_PATH` (a volume path on Railway and the
Macs), migrated by `PRAGMA user_version` (`_MIGRATIONS` maps each version to the
function that builds it, so a new version is one entry plus its test). `MemoryUserStore`
is the desktop default when no path is configured: same behavior, no persistence.

Only metadata lives here, never world bytes -- per-user worlds stay `.bloom` files on
disk (see `docs/save-compatibility.md`), recorded by id in the `worlds` table. The
identity providers (item 2) and Polar webhooks (item 4) are later writers to these
same tables; tier checks read them now. Admin is an override, never stored: pass
`admin=True` for the owner account and local/desktop mode (#560).

A missing user is "free", never an error: route gating (item 3) treats unknown
callers as signed-out viewers. Anything else that names a user who isn't there --
linking, subscribing, or adding a world -- raises `KeyError(user_id)` instead of
inventing one. Recording a world already owned by someone else raises `ValueError`
(it is never silently reassigned); display names are capped at 200 characters and
subscription statuses must be strings. Fresh database files are created owner-only
(`0o600`); a file stamped with a newer schema version than this build reads is
refused instead of opened.
"""

from __future__ import annotations

import os
import sqlite3
import threading
from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Literal

Tier = Literal["free", "paid", "admin"]

# Subscription states that unlock paid. Anything else (past_due, canceled, unpaid,
# incomplete...) is free: #560 downgrades to read/export-only, never silently deletes.
PAID_STATUSES = frozenset({"active", "trialing"})

SCHEMA_VERSION = 2

_ENV_VAR = "MANTLE_BLOOM_USER_DB_PATH"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True)
class User:
    id: int
    created_at: str
    display_name: str | None = None


@dataclass(frozen=True)
class Subscription:
    user_id: int
    status: str
    updated_at: str
    provider_customer_id: str | None = None
    provider_subscription_id: str | None = None
    plan: str = "paid"
    current_period_end: str | None = None


class UserStore(ABC):
    """Users, their linked identities, the subscription cache, and world ownership."""

    @abstractmethod
    def create_user(self, display_name: str | None = None) -> User:
        """Insert a user and return it."""

    @abstractmethod
    def get_user(self, user_id: int) -> User | None:
        """The user, or None for an unknown id (unknown callers are free)."""

    @abstractmethod
    def link_identity(self, provider: str, subject: str, user_id: int, email: str | None = None) -> None:
        """Attach a provider login to a user. Raises `KeyError` for an unknown user and
        `ValueError` if the login is already attached to a different user (merging
        accounts needs an explicit, separately authorized flow -- never a silent
        transfer on link)."""

    @abstractmethod
    def user_for_identity(self, provider: str, subject: str) -> User | None:
        """The user behind a provider login, or None."""

    @abstractmethod
    def set_subscription(
        self,
        user_id: int,
        *,
        status: str,
        provider_customer_id: str | None = None,
        provider_subscription_id: str | None = None,
        plan: str = "paid",
        current_period_end: str | None = None,
    ) -> Subscription:
        """Upsert the cached billing truth for a user. Raises `KeyError` if missing."""

    @abstractmethod
    def get_subscription(self, user_id: int) -> Subscription | None:
        """The cached subscription, or None if the user never subscribed."""

    @abstractmethod
    def add_world(self, world_id: str, user_id: int) -> None:
        """Record ownership of a world file. Raises `KeyError` for an unknown user and
        `ValueError` if the id is already recorded under a different user."""

    @abstractmethod
    def remove_world(self, world_id: str) -> None:
        """Forget a world file. Raises `KeyError` for an unrecorded id."""

    @abstractmethod
    def world_owner(self, world_id: str) -> int | None:
        """The owning user id, or None for an unrecorded world."""

    @abstractmethod
    def user_worlds(self, user_id: int) -> list[str]:
        """Every world id recorded for a user, oldest first (ties broken by world id,
        so both backends agree even when two adds share a timestamp)."""

    @abstractmethod
    def tier(self, user_id: int, *, admin: bool = False) -> Tier:
        """`admin` for the override, `paid` for an active/trialing subscription,
        `free` for everyone else -- including unknown users."""


MAX_DISPLAY_NAME_LEN = 200


def _check_display_name(display_name: str | None) -> None:
    if display_name is not None and len(display_name) > MAX_DISPLAY_NAME_LEN:
        raise ValueError(f"display name longer than {MAX_DISPLAY_NAME_LEN} characters")


def _check_status(status: str) -> None:
    if not isinstance(status, str):
        raise TypeError(f"subscription status must be a string, got {type(status).__name__}")


def _tier_for_status(sub: Subscription | None, admin: bool) -> Tier:
    if admin:
        return "admin"
    if sub is not None and sub.status in PAID_STATUSES:
        return "paid"
    return "free"


class MemoryUserStore(UserStore):
    """In-process store: the desktop default and the unit-test backend."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._users: dict[int, User] = {}
        self._identities: dict[tuple[str, str], tuple[int, str | None]] = {}
        self._subscriptions: dict[int, Subscription] = {}
        self._worlds: dict[str, tuple[int, str]] = {}
        self._next_id = 1

    def create_user(self, display_name: str | None = None) -> User:
        _check_display_name(display_name)
        with self._lock:
            user = User(id=self._next_id, created_at=_now(), display_name=display_name)
            self._users[user.id] = user
            self._next_id += 1
        return user

    def get_user(self, user_id: int) -> User | None:
        return self._users.get(user_id)

    def link_identity(self, provider: str, subject: str, user_id: int, email: str | None = None) -> None:
        with self._lock:
            if user_id not in self._users:
                raise KeyError(user_id)
            hit = self._identities.get((provider, subject))
            if hit is not None and hit[0] != user_id:
                raise ValueError(f"identity {provider}:{subject} is already linked to user {hit[0]}")
            self._identities[(provider, subject)] = (user_id, email)

    def user_for_identity(self, provider: str, subject: str) -> User | None:
        hit = self._identities.get((provider, subject))
        return self._users.get(hit[0]) if hit else None

    def set_subscription(
        self,
        user_id: int,
        *,
        status: str,
        provider_customer_id: str | None = None,
        provider_subscription_id: str | None = None,
        plan: str = "paid",
        current_period_end: str | None = None,
    ) -> Subscription:
        _check_status(status)
        with self._lock:
            if user_id not in self._users:
                raise KeyError(user_id)
            sub = Subscription(
                user_id=user_id,
                status=status,
                updated_at=_now(),
                provider_customer_id=provider_customer_id,
                provider_subscription_id=provider_subscription_id,
                plan=plan,
                current_period_end=current_period_end,
            )
            self._subscriptions[user_id] = sub
        return sub

    def get_subscription(self, user_id: int) -> Subscription | None:
        return self._subscriptions.get(user_id)

    def add_world(self, world_id: str, user_id: int) -> None:
        with self._lock:
            if user_id not in self._users:
                raise KeyError(user_id)
            hit = self._worlds.get(world_id)
            if hit is not None and hit[0] != user_id:
                raise ValueError(f"world {world_id!r} is already owned by user {hit[0]}")
            if hit is None:
                self._worlds[world_id] = (user_id, _now())

    def remove_world(self, world_id: str) -> None:
        with self._lock:
            if world_id not in self._worlds:
                raise KeyError(world_id)
            del self._worlds[world_id]

    def world_owner(self, world_id: str) -> int | None:
        hit = self._worlds.get(world_id)
        return hit[0] if hit else None

    def user_worlds(self, user_id: int) -> list[str]:
        with self._lock:
            owned = [(wid, ts) for wid, (uid, ts) in self._worlds.items() if uid == user_id]
        return [wid for wid, _ in sorted(owned, key=lambda item: (item[1], item[0]))]

    def tier(self, user_id: int, *, admin: bool = False) -> Tier:
        return _tier_for_status(self._subscriptions.get(user_id), admin)


def _migrate_to_1(conn: sqlite3.Connection) -> None:
    """Users, the subscription cache, and world ownership."""
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS users(
            id INTEGER PRIMARY KEY,
            created_at TEXT NOT NULL,
            display_name TEXT
        );
        CREATE TABLE IF NOT EXISTS subscriptions(
            user_id INTEGER PRIMARY KEY REFERENCES users(id),
            provider_customer_id TEXT,
            provider_subscription_id TEXT,
            status TEXT NOT NULL,
            plan TEXT NOT NULL DEFAULT 'paid',
            current_period_end TEXT,
            updated_at TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS worlds(
            world_id TEXT PRIMARY KEY,
            user_id INTEGER NOT NULL REFERENCES users(id),
            created_at TEXT NOT NULL
        );
        CREATE INDEX IF NOT EXISTS idx_worlds_user ON worlds(user_id);
        """
    )


def _migrate_to_2(conn: sqlite3.Connection) -> None:
    """Linked provider logins, for the AuthN item: kept separate so a v1 database
    upgrades cleanly to the current schema."""
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS identities(
            provider TEXT NOT NULL,
            subject TEXT NOT NULL,
            user_id INTEGER NOT NULL REFERENCES users(id),
            email TEXT,
            created_at TEXT NOT NULL,
            PRIMARY KEY (provider, subject)
        );
        CREATE INDEX IF NOT EXISTS idx_identities_user ON identities(user_id);
        """
    )


_MIGRATIONS = {1: _migrate_to_1, 2: _migrate_to_2}


def _migrate(conn: sqlite3.Connection) -> None:
    version = conn.execute("PRAGMA user_version").fetchone()[0]
    if version > SCHEMA_VERSION:
        raise ValueError(
            f"user database version {version} is newer than this build reads ({SCHEMA_VERSION})"
        )
    for target in range(version + 1, SCHEMA_VERSION + 1):
        _MIGRATIONS[target](conn)
        conn.execute(f"PRAGMA user_version = {target}")
    conn.commit()


class SqliteUserStore(UserStore):
    """Durable store on a volume path. One connection plus a lock: every method is
    short, so concurrent requests just queue instead of fighting over the file."""

    def __init__(self, path: str) -> None:
        resolved = os.path.realpath(path)
        fresh = not os.path.exists(resolved)
        parent = os.path.dirname(resolved) or "."
        os.makedirs(parent, exist_ok=True)
        self._conn = sqlite3.connect(resolved, check_same_thread=False)
        self._lock = threading.Lock()
        with self._lock:
            self._conn.execute("PRAGMA foreign_keys = ON")
            _migrate(self._conn)
        if fresh:
            # Emails and provider customer ids live here: keep the file owner-only.
            os.chmod(resolved, 0o600)

    def _one(self, query: str, args: tuple = ()) -> sqlite3.Row | None:
        with self._lock:
            return self._conn.execute(query, args).fetchone()

    def _integrity_owner(self, user_id: int) -> None:
        """Map a foreign-key failure to the missing user, if that is what failed;
        anything else (a NULL status, a broken schema) stays a loud integrity error
        instead of masquerading as an unknown user."""
        with self._lock:
            known = self._conn.execute("SELECT 1 FROM users WHERE id = ?", (user_id,)).fetchone()
        if known is None:
            raise KeyError(user_id)
        raise

    def _run(self, query: str, args: tuple, user_id: int) -> None:
        try:
            with self._lock:
                self._conn.execute(query, args)
                self._conn.commit()
        except sqlite3.IntegrityError:
            self._integrity_owner(user_id)

    def create_user(self, display_name: str | None = None) -> User:
        _check_display_name(display_name)
        with self._lock:
            row = self._conn.execute(
                "INSERT INTO users(display_name, created_at) VALUES (?, ?) RETURNING id, created_at",
                (display_name, _now()),
            ).fetchone()
            self._conn.commit()
        return User(id=row[0], created_at=row[1], display_name=display_name)

    def get_user(self, user_id: int) -> User | None:
        row = self._one("SELECT id, created_at, display_name FROM users WHERE id = ?", (user_id,))
        return User(id=row[0], created_at=row[1], display_name=row[2]) if row else None

    def link_identity(self, provider: str, subject: str, user_id: int, email: str | None = None) -> None:
        # Check and write under one lock hold: the owner check and the upsert must be
        # atomic, and a login attached elsewhere is never silently transferred.
        with self._lock:
            if self._conn.execute("SELECT 1 FROM users WHERE id = ?", (user_id,)).fetchone() is None:
                raise KeyError(user_id)
            row = self._conn.execute(
                "SELECT user_id FROM identities WHERE provider = ? AND subject = ?", (provider, subject)
            ).fetchone()
            if row is not None and row[0] != user_id:
                raise ValueError(f"identity {provider}:{subject} is already linked to user {row[0]}")
            self._conn.execute(
                "INSERT INTO identities(provider, subject, user_id, email, created_at)"
                " VALUES (?, ?, ?, ?, ?)"
                " ON CONFLICT(provider, subject) DO UPDATE SET user_id = excluded.user_id, email = excluded.email",
                (provider, subject, user_id, email, _now()),
            )
            self._conn.commit()

    def user_for_identity(self, provider: str, subject: str) -> User | None:
        row = self._one(
            "SELECT u.id, u.created_at, u.display_name FROM users u"
            " JOIN identities i ON i.user_id = u.id WHERE i.provider = ? AND i.subject = ?",
            (provider, subject),
        )
        return User(id=row[0], created_at=row[1], display_name=row[2]) if row else None

    def set_subscription(
        self,
        user_id: int,
        *,
        status: str,
        provider_customer_id: str | None = None,
        provider_subscription_id: str | None = None,
        plan: str = "paid",
        current_period_end: str | None = None,
    ) -> Subscription:
        _check_status(status)
        sub = Subscription(
            user_id=user_id,
            status=status,
            updated_at=_now(),
            provider_customer_id=provider_customer_id,
            provider_subscription_id=provider_subscription_id,
            plan=plan,
            current_period_end=current_period_end,
        )
        self._run(
            "INSERT INTO subscriptions(user_id, provider_customer_id, provider_subscription_id,"
            " status, plan, current_period_end, updated_at) VALUES (?, ?, ?, ?, ?, ?, ?)"
            " ON CONFLICT(user_id) DO UPDATE SET provider_customer_id = excluded.provider_customer_id,"
            " provider_subscription_id = excluded.provider_subscription_id, status = excluded.status,"
            " plan = excluded.plan, current_period_end = excluded.current_period_end,"
            " updated_at = excluded.updated_at",
            (
                sub.user_id,
                sub.provider_customer_id,
                sub.provider_subscription_id,
                sub.status,
                sub.plan,
                sub.current_period_end,
                sub.updated_at,
            ),
            user_id,
        )
        return sub

    def get_subscription(self, user_id: int) -> Subscription | None:
        row = self._one(
            "SELECT user_id, status, updated_at, provider_customer_id, provider_subscription_id,"
            " plan, current_period_end FROM subscriptions WHERE user_id = ?",
            (user_id,),
        )
        if not row:
            return None
        return Subscription(
            user_id=row[0],
            status=row[1],
            updated_at=row[2],
            provider_customer_id=row[3],
            provider_subscription_id=row[4],
            plan=row[5],
            current_period_end=row[6],
        )

    def add_world(self, world_id: str, user_id: int) -> None:
        # Check and insert under one lock hold: two concurrent first-adds must not
        # both pass the ownership check, and a world recorded under another user is
        # never silently reassigned (that check belongs here, not just the routes).
        with self._lock:
            if self._conn.execute("SELECT 1 FROM users WHERE id = ?", (user_id,)).fetchone() is None:
                raise KeyError(user_id)
            row = self._conn.execute("SELECT user_id FROM worlds WHERE world_id = ?", (world_id,)).fetchone()
            if row is None:
                self._conn.execute(
                    "INSERT INTO worlds(world_id, user_id, created_at) VALUES (?, ?, ?)",
                    (world_id, user_id, _now()),
                )
                self._conn.commit()
            elif row[0] != user_id:
                raise ValueError(f"world {world_id!r} is already owned by user {row[0]}")

    def remove_world(self, world_id: str) -> None:
        with self._lock:
            cur = self._conn.execute("DELETE FROM worlds WHERE world_id = ?", (world_id,))
            self._conn.commit()
            if cur.rowcount == 0:
                raise KeyError(world_id)

    def world_owner(self, world_id: str) -> int | None:
        row = self._one("SELECT user_id FROM worlds WHERE world_id = ?", (world_id,))
        return row[0] if row else None

    def user_worlds(self, user_id: int) -> list[str]:
        with self._lock:
            rows = self._conn.execute(
                "SELECT world_id FROM worlds WHERE user_id = ? ORDER BY created_at, world_id", (user_id,)
            ).fetchall()
        return [row[0] for row in rows]

    def tier(self, user_id: int, *, admin: bool = False) -> Tier:
        return _tier_for_status(self.get_subscription(user_id), admin)


def open_user_store() -> UserStore:
    """The configured store: SQLite at `MANTLE_BLOOM_USER_DB_PATH`, or in-memory
    when unset. A set-but-unusable path fails loudly at startup instead of silently
    running memory-backed and losing every user on restart."""
    path = os.environ.get(_ENV_VAR, "").strip()
    if not path:
        return MemoryUserStore()
    return SqliteUserStore(path)


def seed_test_users(store: UserStore) -> dict[str, int]:
    """One free, one paid and one admin-test user for tier tests. The admin entry is
    used with `admin=True` (admin is an override, never stored); the free user never
    subscribes, the paid one holds an active subscription."""
    free = store.create_user("free-fixture")
    paid = store.create_user("paid-fixture")
    store.set_subscription(
        paid.id, status="active", provider_customer_id="cus_test", provider_subscription_id="sub_test"
    )
    admin = store.create_user("admin-fixture")
    return {"free": free.id, "paid": paid.id, "admin": admin.id}
