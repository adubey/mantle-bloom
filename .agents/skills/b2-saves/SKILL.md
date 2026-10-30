---
name: b2-saves
description: Store, find and move Mantle Bloom save files (.bloom, .mbworld) across machines with a private Backblaze B2 bucket. Use when asked to upload, back up, sync, list, fetch or free up saves, when a save named by an issue or replay isn't on this machine, or to set up B2 on a new machine.
---

# B2 saves

`bin/b2_saves.py` keeps saves in a private B2 bucket. It uses only the Python standard
library, so run it with any `python3` 3.9 or newer, from any checkout.

Recent saves live in both places. The local directory is `MANTLE_BLOOM_SAVES_DIR` (default
`~/Downloads`), the same one [`replay-ab`](../replay-ab/SKILL.md) reads, and saves are not
moved out of it. Bucket copies stay until someone deletes them by hand.

| Task | Command |
|---|---|
| See everything, on any machine | `bin/b2_saves.py list` (`--json` for all metadata) |
| Fetch a save here | `bin/b2_saves.py pull <name>` |
| Upload one save with a note | `bin/b2_saves.py push <file> --note "why it matters"` |
| Upload everything new | `bin/b2_saves.py sync` |
| Preview a local cleanup | `bin/b2_saves.py purge --dry-run` |
| Clean up local copies | `bin/b2_saves.py purge` |
| Protect a local copy | `bin/b2_saves.py pin <name>`, undo with `unpin` |
| Check the setup | `bin/b2_saves.py doctor` |

`list` shows whether each file is `local`, `cloud` or `both`, where `both` means the same
bytes by SHA-1, plus its size, simulated age, upload date, source machine, world name and
note. Metadata is stored with each file and can't be edited later, so give a note at push
time.

## Rules

- **Never put a credential in the repo, an issue, a PR, a log or a chat summary.** The repo
  is public. Credentials stay in `~/.config/mantle-bloom/b2.env`, which is mode 600 and sits
  outside every checkout.
- **Keep the bucket private.** Every push, pull, list and cron run refuses a public bucket.
  `delete-remote` still works, so a mistake can be cleaned up.
- **Pull only from this bucket and treat it as trusted input.** Older `.mbworld` saves are
  pickles, and loading a pickle can run code. SHA-1 checks catch corruption, not tampering.
- **Run `delete-remote` only when the user names that file.** It removes every version from
  the bucket, can't be undone, and needs `--yes`.
- **Don't change the retention settings to get a file deleted.** Ask the user.

## What purge deletes

A local save goes only if all of these hold:

- It matches `mantle-bloom-*.bloom` or `mantle-bloom-*.mbworld`.
- It's older than `B2_SAVES_KEEP_DAYS` (default 14), by file modification time.
- It's outside the `B2_SAVES_KEEP_NEWEST` newest saves (default 5).
- It isn't pinned on this machine.
- A fresh bucket listing shows a file with the same SHA-1.

If the listing fails, nothing is deleted. A pulled save counts as new, so it survives until
it ages out again. Set either `B2_SAVES_KEEP_*` value in `b2.env`.

## Set up a machine

1. Create `~/.config/mantle-bloom/b2.env` with mode 600:

   ```
   B2_APPLICATION_KEY_ID=<key ID>
   B2_APPLICATION_KEY=<secret key, shown once when the key is created>
   B2_BUCKET=<private bucket name>
   ```

   Copy this file between machines over a private channel. Ask the user for values you
   don't have, and don't echo them back. A key limited to this one bucket is better than
   the account's master key. It needs `listBuckets`, `listFiles`, `readFiles` and
   `writeFiles`, plus `deleteFiles` for `delete-remote`.

2. Run `bin/b2_saves.py doctor`, and fix anything it reports.

3. Run `bin/b2_saves.py purge --dry-run`, and show the user what the first run would remove.

4. With the user's go-ahead, run `bin/b2_saves.py install-cron`. It copies the script to
   `~/.local/share/mantle-bloom/`, so the job doesn't depend on a worktree, and runs `sync`
   then `purge` hourly. Run `install-cron` again after the script changes.

5. On macOS, when the saves directory is `~/Downloads`, privacy controls block cron from
   it. Have the user add `/usr/sbin/cron` under System Settings > Privacy & Security > Full
   Disk Access. This grants all cron jobs that access, so mention it first.

6. After the next run, read `~/.local/state/mantle-bloom/b2_saves.log`.
