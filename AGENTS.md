# SFINCS

Before editing any Fortran source under `source/`, read [source/CODING_STYLE.md](source/CODING_STYLE.md) and follow it for every edit: bug fixes, new features and reformatting alike.

Before adding a feature or fixing a bug in the Fortran source, read [source/CODING_DEVELOP.md](source/CODING_DEVELOP.md): how the code is organised, the initialisation and time-loop order, and where new code goes.

[source/CODING_REVIEW.md](source/CODING_REVIEW.md) defines when a Fortran change is complete. Use it in two situations, without waiting to be asked:

- **Before committing or opening a pull request** that touches `source/`: check your changes against every section of the checklist they touch (new variables declared, GPU data, `Makefile.am` order, input formats and defaults, documentation, testing). Fix what you find, and report anything you could not verify or fix.
- **When opening a pull request**: fill in the pull request template as described in section 6 of [source/CODING_DEVELOP.md](source/CODING_DEVELOP.md). Tick items only when they were actually done, tick and add "n/a" when they do not apply, and never tick "I have reviewed my own changes": that box is for the human developer.
- **When reviewing** Fortran changes (a diff, branch or pull request): work through every section the changes touch, in addition to checking the style rules.

## Commits

When an AI coding agent (Claude Code, GitHub Copilot, Codex, etc.) creates a commit, end the commit message with exactly this trailer, after a blank line:

```
Co-Authored-By: SFINCS <sfincs@deltares.nl>
```

Use it instead of the agent's own attribution: do not add `Co-Authored-By` trailers for Claude, Copilot, Codex or any other agent.
