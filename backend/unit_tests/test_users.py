import sqlite3

import pytest
from app import users
from app.users import MemoryUserStore, SqliteUserStore


@pytest.fixture(params=["memory", "sqlite"])
def store(request, tmp_path, monkeypatch):
    if request.param == "memory":
        monkeypatch.delenv("MANTLE_BLOOM_USER_DB_PATH", raising=False)
        return MemoryUserStore()
    return SqliteUserStore(str(tmp_path / "users.db"))


def test_create_and_get_user_round_trips(store):
    user = store.create_user("ada")

    assert store.get_user(user.id) == user
    assert store.get_user(999_999) is None
    assert user.display_name == "ada"
    assert user.created_at


def test_link_and_look_up_identity(store):
    user = store.create_user()
    store.link_identity("github", "123", user.id, email="a@example.com")

    assert store.user_for_identity("github", "123") == user
    assert store.user_for_identity("github", "nope") is None
    assert store.user_for_identity("google", "123") is None


def test_relink_moves_identity_to_new_owner(store):
    first = store.create_user()
    second = store.create_user()
    store.link_identity("github", "123", first.id)
    store.link_identity("github", "123", second.id)

    assert store.user_for_identity("github", "123") == second


def test_identity_for_unknown_user_fails(store):
    with pytest.raises(KeyError):
        store.link_identity("github", "123", 999_999)


def test_subscription_upsert_and_statuses(store):
    ids = users.seed_test_users(store)

    assert store.tier(ids["free"]) == "free"
    assert store.tier(ids["paid"]) == "paid"
    assert store.tier(ids["admin"]) == "free"
    assert store.tier(ids["admin"], admin=True) == "admin"
    assert store.tier(999_999) == "free"

    sub = store.get_subscription(ids["paid"])
    assert sub.status == "active"
    assert sub.provider_customer_id == "cus_test"

    store.set_subscription(ids["paid"], status="trialing")
    assert store.tier(ids["paid"]) == "paid"
    for status in ("past_due", "canceled", "unpaid", "incomplete"):
        store.set_subscription(ids["paid"], status=status)
        assert store.tier(ids["paid"]) == "free"
    assert store.get_subscription(ids["free"]) is None


def test_subscription_for_unknown_user_fails(store):
    with pytest.raises(KeyError):
        store.set_subscription(999_999, status="active")


def test_world_ownership_lifecycle(store):
    ids = users.seed_test_users(store)
    store.add_world("w1", ids["free"])
    store.add_world("w2", ids["free"])
    store.add_world("w3", ids["paid"])

    assert store.world_owner("w1") == ids["free"]
    assert store.world_owner("missing") is None
    assert store.user_worlds(ids["free"]) == ["w1", "w2"]
    assert store.user_worlds(ids["paid"]) == ["w3"]

    store.remove_world("w1")
    assert store.world_owner("w1") is None
    assert store.user_worlds(ids["free"]) == ["w2"]
    with pytest.raises(KeyError):
        store.remove_world("w1")


def test_world_for_unknown_user_fails(store):
    with pytest.raises(KeyError):
        store.add_world("w9", 999_999)


def test_open_user_store_selects_backend(tmp_path, monkeypatch):
    monkeypatch.delenv("MANTLE_BLOOM_USER_DB_PATH", raising=False)
    assert isinstance(users.open_user_store(), MemoryUserStore)

    path = str(tmp_path / "sub" / "users.db")
    monkeypatch.setenv("MANTLE_BLOOM_USER_DB_PATH", path)
    assert isinstance(users.open_user_store(), SqliteUserStore)


def test_open_user_store_fails_loudly_on_unusable_path(tmp_path, monkeypatch):
    blocker = tmp_path / "file"
    blocker.write_text("not a dir")
    monkeypatch.setenv("MANTLE_BLOOM_USER_DB_PATH", str(blocker / "users.db"))

    with pytest.raises(OSError):
        users.open_user_store()


def test_sqlite_schema_version_is_current(tmp_path):
    path = str(tmp_path / "users.db")
    SqliteUserStore(path)

    with sqlite3.connect(path) as conn:
        assert conn.execute("PRAGMA user_version").fetchone()[0] == users.SCHEMA_VERSION


def test_sqlite_v1_upgrades_to_current_with_data(tmp_path):
    path = str(tmp_path / "users.db")
    with sqlite3.connect(path) as conn:
        users._migrate_to_1(conn)
        conn.execute("PRAGMA user_version = 1")
        conn.execute("INSERT INTO users(id, display_name, created_at) VALUES (7, 'old', 't')")
        conn.execute("INSERT INTO worlds(world_id, user_id, created_at) VALUES ('w', 7, 't')")
        conn.commit()

    store = SqliteUserStore(path)

    assert store.get_user(7).display_name == "old"
    assert store.world_owner("w") == 7
    store.link_identity("github", "1", 7)
    assert store.user_for_identity("github", "1").id == 7
    with sqlite3.connect(path) as conn:
        assert conn.execute("PRAGMA user_version").fetchone()[0] == users.SCHEMA_VERSION


def test_sqlite_reopen_keeps_data_and_version(tmp_path):
    path = str(tmp_path / "users.db")
    first = SqliteUserStore(path)
    ids = users.seed_test_users(first)
    first.add_world("w", ids["paid"])

    second = SqliteUserStore(path)

    assert second.get_user(ids["paid"]) == first.get_user(ids["paid"])
    assert second.tier(ids["paid"]) == "paid"
    assert second.world_owner("w") == ids["paid"]


def test_unknown_tier_user_is_free_without_backend():
    assert MemoryUserStore().tier(123) == "free"
    assert MemoryUserStore().tier(123, admin=True) == "admin"


def test_seed_test_users_shape():
    ids = users.seed_test_users(MemoryUserStore())

    assert set(ids) == {"free", "paid", "admin"}
    assert len(set(ids.values())) == 3
