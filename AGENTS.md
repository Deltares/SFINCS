# SFINCS

Before editing any Fortran source under `source/`, read [source/CODING_STYLE.md](source/CODING_STYLE.md) and follow it for every edit: bug fixes, new features and reformatting alike.

## Commits

When an AI coding agent (Claude Code, GitHub Copilot, Codex, etc.) creates a commit, end the commit message with exactly this trailer, after a blank line:

```
Co-Authored-By: SFINCS <sfincs@deltares.nl>
```

Use it instead of the agent's own attribution: do not add `Co-Authored-By` trailers for Claude, Copilot, Codex or any other agent.
