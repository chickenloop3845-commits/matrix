# Agent rules (LLM)

This file contains rules intended for LLM/code-agent sessions. If you are an assistant working in this repo, **read this file first**.

## Git workflow rules

1) **Before committing/pushing**

- Summarize what you are about to commit (key files touched + intent).
- Verify the change set does **not** appear to contain secrets (tokens, private keys, credentials, etc.). If it does, do not commit.

2) **Before pushing**

- Confirm you are **not** in a detached HEAD state.
  - `git status -sb` must **not** show `HEAD (no branch)`.
  - If detached, fix it **before** pushing (e.g., switch to the intended branch and cherry-pick the detached commit(s)).

3) **After pushing**

- If any PDFs or Markdown files were added/updated, provide the **direct GitHub link(s)** to those files (e.g., `https://github.com/<org>/<repo>/blob/<branch>/<path>.pdf` or `https://github.com/<org>/<repo>/blob/<branch>/<path>.md`).
