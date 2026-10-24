# Run mantle-bloom locally

This guide covers running the app from a checkout and enabling GitHub sign-in on localhost.
Local runs skip authentication by default; configure OAuth and turn off that bypass when you
want to test the real sign-in flow.

## Configure GitHub OAuth

Create a separate OAuth App for local development under **GitHub → Settings → Developer
settings → OAuth Apps → New OAuth App**. GitHub's [OAuth app creation guide](https://docs.github.com/en/apps/oauth-apps/building-oauth-apps/creating-an-oauth-app)
describes the registration fields.

For a server running on port `8425`, enter:

| OAuth App field | Value |
|---|---|
| Homepage URL | `http://localhost:8425` |
| Authorization callback URL | `http://localhost:8425/auth/github/callback` |

Keep the port in this callback URL the same as the port used to start the server. GitHub
redirects to the app's authorization callback after sign-in; see its [callback URL guide](https://docs.github.com/en/apps/creating-github-apps/registering-a-github-app/about-the-user-authorization-callback-url).

After creating the OAuth App, copy its **Client ID** and generate a **Client secret**. The
server also needs the numeric GitHub ID of the one account allowed to use the full app. Get
the ID for the currently authenticated `gh` account with `gh api user -q .id`.

Set the values in the shell that will run the server:

```bash
export MANTLE_BLOOM_GITHUB_CLIENT_ID='<client ID>'
export MANTLE_BLOOM_GITHUB_CLIENT_SECRET='<client secret>'
export MANTLE_BLOOM_GITHUB_OWNER_ID="$(gh api user -q .id)"
export MANTLE_BLOOM_SESSION_SECRET="$(python -c 'import secrets; print(secrets.token_urlsafe(48))')"
export MANTLE_BLOOM_LOCAL_NO_LOGIN=0
```

The session secret must be at least 32 bytes. Keep both secrets out of source control. With
`MANTLE_BLOOM_LOCAL_NO_LOGIN=0`, localhost uses the same sign-in check as a remote request;
without it, the default local bypass treats localhost as the owner.

## Start the app

Install dependencies once, then start the single-process app on the same port used in the
OAuth callback:

```bash
./bin/setup.sh
./bin/restart.sh --port 8425
```

Open [http://localhost:8425/login](http://localhost:8425/login) to sign in. The full app is at
[http://localhost:8425/](http://localhost:8425/). Generate or load a world there, then preview
the logged-out viewer at [http://localhost:8425/demo/view](http://localhost:8425/demo/view).
The demo viewer follows the guest step allowance and advances the world automatically.

`/demo/` opens the login page even when the local no-login bypass is enabled. Use it to preview
the login screen; to complete GitHub sign-in, configure the OAuth values above.

If port `8425` is already in use, choose another available port and update both the OAuth App's
homepage and callback URLs to match.
