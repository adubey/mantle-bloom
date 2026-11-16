# Issue 920: penetration test of a local copy

Scope, checks and results for the penetration-test part of
[#920](https://github.com/adubey/mantle-bloom/issues/920). The dependency scan and frontend
XSS audit are in this folder's [README](README.md).

## Scope

- **Target:** the application code, exercised in-process and on a loopback copy of the
  hosted entry point with throwaway credentials. The owner took the deployed service out of
  scope, so nothing was sent to it.
- **Method:** each area from #920 is checked by unit tests that stay in the suite, so the
  result is rechecked on every pull request. Tests that existed before this pass are
  credited where they already answer the question.
- **Not covered:** anything only the deployment can show (TLS and edge behaviour, cookie
  flags as served, whether its variables match the repository), sustained load, and the
  base-image scan ([#926](https://github.com/adubey/mantle-bloom/issues/926)).

## Results

One new finding (S12, fixed with these tests). S1, unbounded request bodies, was found
earlier in this pass and fixed in PR #924.

| Area | Check | Result | Test |
|---|---|---|---|
| OAuth | Start sends `state` and an S256 PKCE challenge; callback needs the matching signed state cookie; a missing or mismatched state fails | Holds | `test_auth.py`: `test_login_start_redirects_to_github_with_state_and_pkce`, `test_callback_failures_redirect_to_login`, `test_callback_without_state_cookie_fails` |
| OAuth | A finished sign-in replayed in the same browser | Fails: the callback clears the state cookie. A copy of that cookie taken beforehand stays valid for its 10 minutes (it is stateless, like the session), but is useless without a fresh GitHub code for the same PKCE verifier | new: `test_a_finished_sign_in_cannot_be_replayed_in_the_same_browser` |
| OAuth | Another GitHub account gets no session | Holds | `test_other_github_account_gets_no_session` |
| Session forgery | Wrong secret, changed owner, expired token | Refused | `test_forged_or_stale_sessions_are_rejected`, `test_expired_token_is_rejected` |
| Session forgery | A guest-allowance or OAuth-state token, signed with the same secret and naming the owner, used as a session | Refused | new: `test_tokens_minted_for_another_purpose_are_not_sessions` |
| Session forgery | Malformed cookies, including correctly signed bodies that aren't a well-formed session (no expiry, wrong types) | Refused with `401`, never a `500` | new: `test_malformed_session_cookies_are_refused_not_errors` |
| Session replay | A copied cookie after logout | Still valid until expiry or secret rotation: accepted (S10) and now pinned | new: `test_copied_session_outlives_logout_until_it_expires` |
| CSRF | Cross-site `Origin`, `Origin: null`, `Sec-Fetch-Site: cross-site`; the dev origin against a remote server | Refused | `test_cross_site_writes_are_refused`, `test_dev_frontend_origin_may_not_write_remotely`, `test_public_view.py::test_cross_site_step_is_refused` |
| Local bypass | Loopback request relayed by a proxy; hosted entry point on loopback | No admin | `test_proxied_loopback_requests_need_login`, `test_hosted.py::test_auth_guard_runs_before_proxy_and_worker_routes_are_private` |
| Authorization | Every route of the simulator app is owner-only unless on an explicit public list; every owner route answers `401` to every method without a session | Holds | new: `test_auth.py::test_every_route_is_owner_only_or_listed_as_public` |
| Authorization | The same inventory for the hosted public process, without waking the worker | Holds | new: `test_hosted.py::test_every_hosted_route_is_owner_only_or_listed_as_public` |
| Authorization | A root-path prefix can't slip past the guard | Holds | `test_prefixed_app_still_guards_world_routes` |
| Authorization | Private worlds on `/public/*` | `404` | `test_public_view.py` (several) |
| Worker proxy | `/_worker/*` from outside | `404` | `test_auth_guard_runs_before_proxy_and_worker_routes_are_private` |
| Worker proxy | The worker routes the same path the public process matched | **Did not hold (S12), fixed** | new: `test_proxy_forwards_the_path_this_process_matched`, `test_proxy_refuses_a_path_the_client_library_would_rewrite`, `test_refused_path_releases_the_worker_lease` |
| Worker proxy | The proxy keeps no cookies between callers | Holds | `test_proxy_does_not_keep_cookies_or_hide_a_snapshot_when_the_worker_exits` |
| Guest allowance | Limit enforced; forged cookie ignored; stale or shared steps not charged | Holds. Clearing the cookie resets it: accepted (#555 F4) | `test_public_view.py`: `test_allowance_runs_out_after_the_limit`, `test_forged_allowance_cookie_is_ignored` |
| Resource bounds | Render dimensions, step years, animation step years, run step years | Refused past the bound | `test_main.py`: `test_render_out_of_range_dimensions_return_400`, `test_step_rejects_invalid_years`, `test_animate_rejects_invalid_step_years`, `test_run_rejects_step_years_over_the_max` |
| Resource bounds | Request body size | Did not hold (S1), fixed in PR #924 | `test_body_limit.py` |
| Resource bounds | Concurrent worlds | Capped | `test_world_registry.py` |
| Upload | Over the 512 MB cap; a legacy pickle | `413`; refused before unpickling | `test_main.py`: `test_load_rejects_an_upload_over_the_size_cap`, `test_uploads_refuse_a_legacy_pickled_save_without_unpickling_it` |
| Upload | A small save that expands without limit; other malformed v5 input | Open: [#925](https://github.com/adubey/mantle-bloom/issues/925) | none yet |

### S12: the worker could route a different path from the one authorised (fixed)

Severity: Low. Surface: `hosted_worker._upstream_request`.

The public process forwarded the raw request path to the worker through `httpx`, which
removes `.` and `..` segments. The public process itself does not, so the two could
disagree. Measured in-process: a request the public process matched as
`/public/worlds/{world_id}/step` with a world id of `..` was sent to the worker as
`/public/step`.

No harmful case was found. A world id can't contain a slash, so a logged-out visitor
could move the path up by one segment only, to worker paths that don't exist. The owner's
catch-all route could reach the worker's private `/_worker/*` routes this way, but the
owner controls the worker anyway. It is fixed because the guard's guarantee is that the
worker routes what the public process authorised, and that should not depend on which
routes happen to exist.

Fix: the proxy answers `404` when the path it would send differs from the path it
received. Percent-encoded paths are forwarded unchanged, and both processes decode them the
same way.
