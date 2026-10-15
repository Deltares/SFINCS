# Line endings in the SFINCS repository

All text files in this repository are checked out with Unix (LF) line endings on every platform.
This is set by the `.gitattributes` file in the repository root and overrides your personal git
settings (such as `core.autocrlf` on Windows).

The reason: shell scripts, `Makefile.am` and other autotools files fail with CRLF line endings,
for example when a Windows checkout is built in WSL or Docker.

Windows-only files keep Windows (CRLF) line endings: `*.bat`, `*.sln`, `*.vfproj`, `*.vcxproj`.

## Updating an existing clone

A fresh clone needs nothing. An existing clone does not convert by pulling alone, because git only
rewrites a file on disk when its content changes. Convert it once by hand:

```
git status                  # must be clean: commit or stash your work first
git rm --cached -r -q .
git reset --hard
```

`git reset --hard` discards uncommitted changes to tracked files, so make sure the working tree is
clean first. Untracked files (build output, model runs) are left alone.

## Checking the result

```
git ls-files --eol | grep "w/crlf"
```

This should list only the `.bat`, `.sln`, `.vfproj` and `.vcxproj` files.

## Good to know

- Branches created before `.gitattributes` was added do not have this rule. Switching to one can
  leave you with a mix of line endings until `main` is merged into that branch.
- Set your editor to save new files with LF. Git converts them to LF on commit either way, but it
  avoids warnings and a mixed working copy.
- If a new Windows-only file type needs CRLF, add it to `.gitattributes`.
