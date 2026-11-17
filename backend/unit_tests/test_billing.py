import base64
import hashlib
import hmac
import json
import time
import urllib.parse
from types import SimpleNamespace

import pytest
from fastapi.testclient import TestClient

from app import auth, billing, users
from app.main import app

SECRET = "s" * 40
RAW_WEBHOOK_SECRET = "polar_whs_test_secret_1234567890"
B64_WEBHOOK_SECRET = base64.b64encode(RAW_WEBHOOK_SECRET.encode()).decode()
SAME_SITE = {"Origin": "https://testserver"}


@pytest.fixture(autouse=True)
def configured(monkeypatch):
    monkeypatch.setenv("MANTLE_BLOOM_SESSION_SECRET", SECRET)
    monkeypatch.setenv("MANTLE_BLOOM_GOOGLE_CLIENT_ID", "google-id")
    monkeypatch.setenv("MANTLE_BLOOM_GOOGLE_CLIENT_SECRET", "google-secret")
    monkeypatch.setenv("POLAR_ACCESS_TOKEN", "polar_oat_test")
    monkeypatch.setenv("POLAR_WEBHOOK_SECRET", B64_WEBHOOK_SECRET)
    monkeypatch.setenv("POLAR_PRODUCT_ID", "prod_123")
    monkeypatch.setenv("POLAR_ENVIRONMENT", "sandbox")
    monkeypatch.delenv("MANTLE_BLOOM_PUBLIC_URL", raising=False)
    monkeypatch.delenv("MANTLE_BLOOM_LOCAL_NO_LOGIN", raising=False)
    monkeypatch.setattr(users, "_shared", users.MemoryUserStore())


@pytest.fixture
def remote():
    return TestClient(app, base_url="https://testserver", follow_redirects=False)


def _sign_in_google(client: TestClient, monkeypatch, sub: str = "1077", name: str = "Ada Lovelace") -> users.User:
    def fake(req):
        if req.full_url == auth.GOOGLE_TOKEN_URL:
            return {"access_token": "ya29.test"}
        return {"sub": sub, "name": name, "email": "ada@example.com", "email_verified": True}

    monkeypatch.setattr(auth, "_http_json", fake)
    start = client.get("/auth/google")
    state = urllib.parse.parse_qs(urllib.parse.urlsplit(start.headers["location"]).query)["state"][0]
    resp = client.get("/auth/google/callback", params={"code": "abc", "state": state})
    assert resp.status_code == 303
    user = users.shared_store().user_for_identity("google", sub)
    assert user is not None
    return user


def _signed_webhook(payload: dict, secret: bytes = RAW_WEBHOOK_SECRET.encode(), *, msg_id: str = "msg_1") -> tuple[bytes, dict]:
    body = json.dumps(payload).encode()
    ts = str(int(time.time()))
    signed_content = f"{msg_id}.{ts}.{body.decode()}".encode()
    sig = base64.b64encode(hmac.new(secret, signed_content, hashlib.sha256).digest()).decode()
    headers = {
        "Content-Type": "application/json",
        "webhook-id": msg_id,
        "webhook-timestamp": ts,
        "webhook-signature": f"v1,{sig}",
    }
    return body, headers


def _subscription_event(
    event_type: str,
    *,
    status: str,
    external_id: str | None,
    customer_id: str = "cus_1",
    subscription_id: str = "sub_1",
    current_period_end: str | None = "2026-12-01T00:00:00Z",
) -> dict:
    return {
        "type": event_type,
        "timestamp": "2026-11-01T00:00:00Z",
        "api_version": "2026-10",
        "data": {
            "id": subscription_id,
            "created_at": "2026-11-01T00:00:00Z",
            "modified_at": None,
            "amount": 500,
            "currency": "usd",
            "recurring_interval": "month",
            "recurring_interval_count": 1,
            "status": status,
            "current_period_start": "2026-11-01T00:00:00Z",
            "current_period_end": current_period_end,
            "current_meter_period_start": None,
            "current_meter_period_end": None,
            "trial_start": None,
            "trial_end": None,
            "cancel_at_period_end": False,
            "canceled_at": None,
            "started_at": "2026-11-01T00:00:00Z",
            "ends_at": None,
            "ended_at": None,
            "past_due_at": None,
            "pause_at_period_end": False,
            "paused_at": None,
            "resumes_at": None,
            "customer_id": customer_id,
            "product_id": "prod_123",
            "discount_id": None,
            "checkout_id": "chk_1",
            "seats": None,
            "units": None,
            "customer_cancellation_reason": None,
            "customer_cancellation_comment": None,
            "metadata": {},
            "custom_field_data": {},
            "customer": {
                "id": customer_id,
                "created_at": "2026-11-01T00:00:00Z",
                "modified_at": None,
                "metadata": {},
                "external_id": external_id,
                "email": "ada@example.com",
                "email_verified": True,
                "type": "individual",
                "name": "Ada",
                "billing_name": None,
                "billing_address": None,
                "tax_id": None,
                "locale": None,
                "organization_id": "org_1",
                "default_payment_method_id": None,
                "deleted_at": None,
                "first_user_event_at": None,
                "avatar_url": None,
            },
            "product": {
                "id": "prod_123",
                "created_at": "2026-11-01T00:00:00Z",
                "modified_at": None,
                "trial_interval": None,
                "trial_interval_count": None,
                "name": "Mantle-Bloom Pro",
                "description": None,
                "visibility": "public",
                "recurring_interval": "month",
                "recurring_interval_count": 1,
                "is_recurring": True,
                "is_archived": False,
                "is_deletable": False,
                "organization_id": "org_1",
                "metadata": {},
                "prices": [],
                "benefits": [],
                "medias": [],
                "attached_custom_fields": [],
                "meter_interval": None,
                "meter_interval_count": None,
            },
            "discount": None,
            "prices": [],
            "meters": [],
            "pending_update": None,
        },
    }


def test_checkout_requires_signed_in_user(remote):
    resp = remote.post("/billing/checkout", headers=SAME_SITE)
    assert resp.status_code == 401


def test_portal_requires_signed_in_user(remote):
    resp = remote.post("/billing/portal", headers=SAME_SITE)
    assert resp.status_code == 401


def test_checkout_creates_session_and_redirects(remote, monkeypatch):
    user = _sign_in_google(remote, monkeypatch)
    calls = []

    class FakeClient:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        class checkouts:
            @staticmethod
            def create(**kwargs):
                calls.append(kwargs)
                return SimpleNamespace(url="https://sandbox.polar.sh/checkout/chk_1")

    monkeypatch.setattr(billing, "_client", lambda config: FakeClient())
    resp = remote.post("/billing/checkout", headers=SAME_SITE)
    assert resp.status_code == 303
    assert resp.headers["location"] == "https://sandbox.polar.sh/checkout/chk_1"
    assert calls == [{
        "products": ["prod_123"],
        "external_customer_id": str(user.id),
        "success_url": "https://testserver/view?billing=success",
        "return_url": "https://testserver/view",
    }]


def test_portal_creates_session_and_redirects(remote, monkeypatch):
    user = _sign_in_google(remote, monkeypatch)
    calls = []

    class FakeClient:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        class customer_sessions:
            @staticmethod
            def create(**kwargs):
                calls.append(kwargs)
                return SimpleNamespace(customer_portal_url="https://sandbox.polar.sh/portal/sess_1")

    monkeypatch.setattr(billing, "_client", lambda config: FakeClient())
    resp = remote.post("/billing/portal", headers=SAME_SITE)
    assert resp.status_code == 303
    assert resp.headers["location"] == "https://sandbox.polar.sh/portal/sess_1"
    assert calls == [{
        "external_customer_id": str(user.id),
        "return_url": "https://testserver/view",
    }]


def test_cross_site_billing_posts_are_refused(remote, monkeypatch):
    _sign_in_google(remote, monkeypatch)
    assert remote.post("/billing/checkout", headers={"Origin": "https://evil.example"}).status_code == 403
    assert remote.post("/billing/portal", headers={"Origin": "https://evil.example"}).status_code == 403


def test_webhook_rejects_missing_or_invalid_signature(remote):
    body = json.dumps(_subscription_event("subscription.active", status="active", external_id="1")).encode()
    assert remote.post("/billing/webhook", content=body).status_code == 403

    _, bad_headers = _subscription_headers_with_wrong_secret(body)
    assert remote.post("/billing/webhook", content=body, headers=bad_headers).status_code == 403


def _subscription_headers_with_wrong_secret(body: bytes) -> tuple[bytes, dict]:
    return _signed_webhook(json.loads(body), secret=b"wrong-secret")


def test_webhook_lifecycle_updates_user_store_and_preserves_records(remote, monkeypatch):
    user = _sign_in_google(remote, monkeypatch)
    store = users.shared_store()
    store.add_world("world-1", user.id)
    assert store.tier(user.id) == "free"

    # Activation flips tier to paid and exposes it in /auth/session.
    body, headers = _signed_webhook(_subscription_event("subscription.active", status="active", external_id=str(user.id)))
    assert remote.post("/billing/webhook", content=body, headers=headers).status_code == 202
    assert store.tier(user.id) == "paid"
    sub = store.get_subscription(user.id)
    assert sub is not None
    assert sub.provider_customer_id == "cus_1"
    assert sub.provider_subscription_id == "sub_1"
    assert sub.current_period_end == "2026-12-01T00:00:00Z"
    assert remote.get("/auth/session").json()["user"]["tier"] == "paid"

    # Past-due downgrades to free without deleting user or world records; admin override unaffected.
    body, headers = _signed_webhook(_subscription_event("subscription.past_due", status="past_due", external_id=str(user.id)))
    assert remote.post("/billing/webhook", content=body, headers=headers).status_code == 202
    assert store.tier(user.id) == "free"
    assert store.tier(user.id, admin=True) == "admin"
    assert store.user_worlds(user.id) == ["world-1"]

    # Re-activation and revocation.
    body, headers = _signed_webhook(_subscription_event("subscription.updated", status="active", external_id=str(user.id)))
    assert remote.post("/billing/webhook", content=body, headers=headers).status_code == 202
    assert store.tier(user.id) == "paid"

    body, headers = _signed_webhook(_subscription_event("subscription.revoked", status="canceled", external_id=str(user.id)))
    assert remote.post("/billing/webhook", content=body, headers=headers).status_code == 202
    assert store.tier(user.id) == "free"
    assert store.get_user(user.id) is not None
    assert store.user_worlds(user.id) == ["world-1"]


def test_subscription_survives_sqlite_restart(tmp_path, monkeypatch):
    db_path = tmp_path / "users.db"
    monkeypatch.setattr(users, "_shared", users.SqliteUserStore(db_path))
    remote = TestClient(app, base_url="https://testserver", follow_redirects=False)
    user = _sign_in_google(remote, monkeypatch)

    body, headers = _signed_webhook(_subscription_event("subscription.active", status="active", external_id=str(user.id)))
    assert remote.post("/billing/webhook", content=body, headers=headers).status_code == 202

    # Simulate service restart by reopening SqliteUserStore from the same path.
    reopened = users.SqliteUserStore(db_path)
    assert reopened.tier(user.id) == "paid"
    sub = reopened.get_subscription(user.id)
    assert sub is not None and sub.provider_subscription_id == "sub_1"


def test_unconfigured_billing_fails_closed_while_local_stays_unlocked(remote, monkeypatch):
    _sign_in_google(remote, monkeypatch)
    monkeypatch.delenv("POLAR_ACCESS_TOKEN")
    monkeypatch.delenv("POLAR_WEBHOOK_SECRET")
    monkeypatch.delenv("POLAR_PRODUCT_ID")

    assert remote.get("/auth/session").json()["billing_available"] is False
    assert remote.post("/billing/checkout", headers=SAME_SITE).status_code == 503
    assert remote.post("/billing/portal", headers=SAME_SITE).status_code == 503
    assert remote.post("/billing/webhook", content=b"{}").status_code == 503

    local = TestClient(app, base_url="http://localhost", client=("127.0.0.1", 50000), follow_redirects=False)
    session = local.get("/auth/session").json()
    assert session["admin"] is True
    assert session["local"] is True
    assert local.get("/worlds").status_code == 200
