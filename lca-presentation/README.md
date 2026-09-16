## LCA Presentation

Quarto/revealjs slide deck on Lowest Common Ancestor, with an
application in phylogenetic trees.

### Build

```
quarto render main.qmd
```

Or preview live while editing:

```
quarto preview main.qmd
```

### Structure

- `main.qmd` — top-level deck, includes every section in order.
  Owned by the leader.
- `slides/` — one file per section, one owner per file.
- `code/lca.py` — shared LCA implementation, imported by
  `slides/03-python-impl.qmd`.
- `assets/images/` — diagrams and figures.
- `assets/data/` — any sample data (e.g. tree definitions) shared
  across sections.

### Task ownership (initial block — subject to leader review)

| File | Section | Owner |
|---|---|---|
| `slides/01-problem-idea.qmd` | Problem definition + core idea | Teammate 1 |
| `slides/02-hand-trace.qmd` | Hand-traced example | Teammate 2 |
| `slides/03-python-impl.qmd` | Python implementation | Teammate 3 |
| `slides/04-visualization.qmd` | Visualization | Teammate 4 |
| `slides/05-complexity.qmd` | Time + space complexity | Teammate 5 |
| `slides/06-application.qmd` | Application (phylogenetics + LeetCode) | Leader |

### Key dependency

Teammate 3's choice of LCA variant (naive / binary lifting / Euler
tour + sparse table) drives sections 2, 4, and 5. Decide this first
and communicate it to the team before everyone starts writing.

### Rendering individual slides

Each teammate can render just their own section while working on
it, without needing `main.qmd` to be finished:

```
quarto render slides/03-python-impl.qmd
```

This works because `_quarto.yml` at the root makes this a Quarto
*project*, and the project's `format: revealjs` setting applies
even to a lone file rendered this way — no per-file YAML needed.

### Gotcha: paths in included files

Because `main.qmd` pulls in every file under `slides/` via
`{{< include >}}`, code that references `code/` or
`assets/images/` needs a path that works **both** when the file is
rendered on its own **and** when it's included into `main.qmd`
(which changes the execution working directory to the project
root). Plain relative paths like `../code` or `code` will each
break in one of the two modes.

The fix used throughout this project: Quarto sets the
`QUARTO_PROJECT_DIR` environment variable to the project's
absolute path on every render, regardless of which file triggered
it. Build paths off that instead of a relative path:

```python
import os
project_dir = os.environ.get("QUARTO_PROJECT_DIR", os.getcwd())
path_to_code = os.path.join(project_dir, "code")
```

See `slides/03-python-impl.qmd` and `slides/04-visualization.qmd`
for working examples. If you add new files that read from `code/`
or `assets/`, follow the same pattern rather than a relative path.

### What not to commit

`.gitignore` excludes rendered output (`*.html`, `*_files/`,
`.quarto/`) and Python/OS junk. Only source files — `.qmd`, `.py`,
images/data you actually authored — should be committed. Rendered
output is regenerated locally by whoever runs `quarto render`, so
committing it just causes merge conflicts and repo bloat without
adding anything reproducible.

If you later want a hosted copy (e.g. GitHub Pages) for the final
version, that's a separate step — publish the rendered output to a
`gh-pages` branch or a CI job, not the main branch.

### Git workflow

- Everyone commits to their own file under `slides/` — no merge
  conflicts on content.
- Only touch `main.qmd` or `_quarto.yml` if you're changing the
  overall structure — coordinate with the leader first.
