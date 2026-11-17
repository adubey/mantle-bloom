"""Polar subscription billing: checkout, customer portal, and webhook sync (#945).

Remote billing routes fail closed (503) when Polar credentials are unconfigured.
Local/desktop loopback mode stays unlocked without any Polar configuration.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Literal

from fastapi import APIRouter, HTTPException, Request, Response
from fastapi.responses import RedirectResponse
from polar import PolarError
from polar.v2026_10 import Polar
from polar.v2026_10.webhooks import (
    PolarWebhookError,
    PolarWebhookUnknownTypeError,
    PolarWebhookVerificationError,
    validate_event,
)

from . import auth, users

_logger = logging.getLogger(__name__)
_VALID_ENVIRONMENTS = frozenset({"production", "sandbox"})

router = APIRouter()


@dataclass(frozen=True)
class BillingConfig:
    access_token: str
    webhook_secret: str
    product_id: str
    environment: Literal["sandbox", "production"] = "production"
    environment_valid: bool = True

    @property
    def checkout_configured(self) -> bool:
        return bool(self.access_token and self.product_id and self.environment_valid)

    @property
    def webhook_configured(self) -> bool:
        return bool(self.webhook_secret)

    @property
    def configured(self) -> bool:
        return self.checkout_configured and self.webhook_configured


def load_config() -> BillingConfig:
    """Read Polar billing settings from the environment on every call."""
    raw_env = os.environ.get("POLAR_ENVIRONMENT", "").strip().lower()
    env_valid = not raw_env or raw_env in _VALID_ENVIRONMENTS
    environment: Literal["sandbox", "production"] = "sandbox" if raw_env == "sandbox" else "production"
    return BillingConfig(
        access_token=os.environ.get("POLAR_ACCESS_TOKEN", "").strip(),
        webhook_secret=os.environ.get("POLAR_WEBHOOK_SECRET", "").strip(),
        product_id=os.environ.get("POLAR_PRODUCT_ID", "").strip(),
        environment=environment,
        environment_valid=env_valid,
    )


def _client(config: BillingConfig) -> Polar:
    return Polar(config.access_token, environment=config.environment)


def _base_url(request: Request, auth_config: auth.AuthConfig) -> str:
    return auth_config.public_url or str(request.base_url).rstrip("/")


def _require_signed_in(request: Request) -> tuple[users.User, auth.AuthConfig]:
    auth_config = auth.load_config(request)
    signed_in = auth.signed_in_user(request, auth_config)
    if signed_in is None:
        raise HTTPException(401, "Sign-in required")
    return signed_in[0], auth_config


@router.post("/billing/checkout", include_in_schema=False)
def create_checkout(request: Request) -> RedirectResponse:
    user, auth_config = _require_signed_in(request)
    config = load_config()
    if not config.checkout_configured:
        raise HTTPException(503, "Billing is not configured")
    base = _base_url(request, auth_config)
    try:
        with _client(config) as client:
            checkout = client.checkouts.create(
                products=[config.product_id],
                external_customer_id=str(user.id),
                success_url=f"{base}/view?billing=success",
                return_url=f"{base}/view",
            )
    except PolarError as exc:
        _logger.warning("Polar checkout creation failed: %s", exc)
        raise HTTPException(502, "Could not start checkout session") from exc
    if not checkout.url:
        raise HTTPException(502, "Could not start checkout session")
    return RedirectResponse(checkout.url, status_code=303)


@router.post("/billing/portal", include_in_schema=False)
def create_portal(request: Request) -> RedirectResponse:
    user, auth_config = _require_signed_in(request)
    config = load_config()
    if not config.checkout_configured:
        raise HTTPException(503, "Billing is not configured")
    base = _base_url(request, auth_config)
    try:
        with _client(config) as client:
            session = client.customer_sessions.create(
                external_customer_id=str(user.id),
                return_url=f"{base}/view",
            )
    except PolarError as exc:
        _logger.warning("Polar customer portal creation failed: %s", exc)
        raise HTTPException(502, "Could not open billing portal") from exc
    if not session.customer_portal_url:
        raise HTTPException(502, "Could not open billing portal")
    return RedirectResponse(session.customer_portal_url, status_code=303)


def _apply_subscription_event(event_type: str, sub: object) -> None:
    customer = getattr(sub, "customer", None)
    external_id = getattr(customer, "external_id", None) if customer is not None else None
    if not isinstance(external_id, str) or not external_id.isdigit():
        _logger.warning("Ignoring Polar %s without numeric external_id", event_type)
        return
    user_id = int(external_id)
    raw_status = str(getattr(sub, "status", ""))
    status = "canceled" if event_type == "subscription.revoked" and raw_status in users.PAID_STATUSES else raw_status
    try:
        users.shared_store().set_subscription(
            user_id,
            status=status,
            provider_customer_id=getattr(sub, "customer_id", None),
            provider_subscription_id=getattr(sub, "id", None),
            plan="paid",
            current_period_end=getattr(sub, "current_period_end", None),
        )
    except KeyError:
        _logger.warning("Ignoring Polar %s for unknown user %s", event_type, user_id)


@router.post("/billing/webhook", include_in_schema=False)
async def polar_webhook(request: Request) -> Response:
    config = load_config()
    if not config.webhook_configured:
        raise HTTPException(503, "Billing is not configured")
    body = await request.body()
    try:
        event = validate_event(body, dict(request.headers), config.webhook_secret)
    except PolarWebhookVerificationError as exc:
        raise HTTPException(403, "Invalid webhook signature") from exc
    except PolarWebhookUnknownTypeError:
        return Response(status_code=202)
    except PolarWebhookError as exc:
        raise HTTPException(400, "Invalid webhook payload") from exc

    if event.type.startswith("subscription."):
        _apply_subscription_event(event.type, event.data)
    return Response(status_code=202)
