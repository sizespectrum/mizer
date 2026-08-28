# Building documentation and the mizer website

Use this skill when updating package documentation, adding or renaming exported
functions, editing agent skills / guide vignettes, or generating website
artifacts.

## Auto-generated files — never edit directly

- `NAMESPACE`, `man/*.Rd` — `devtools::document()`
- `R/RcppExports.R`, `src/RcppExports.cpp` — `devtools::document()` / `Rcpp::compileAttributes()`
- `vignettes/guide-*.Rmd` — from `inst/skills/*/SKILL.md` via `dev_scripts/build_guides.R`
- `vignettes/upgrading.Rmd` — from `inst/skills/upgrade-mizer-code/SKILL.md` plus its `references/mizer-*.md` files, appended newest release first
- `docs/llms.txt`, `inst/llms.txt` — `dev_scripts/build_llms.R`
- `docs/`, and any `vignettes/*.html` — artefacts of the website build (section 3)

---

## 1. Guide vignettes and agent skills

Agent skills in `inst/skills/<topic>/SKILL.md` are the single source of truth for both AI agents and the pkgdown guide articles in `vignettes/guide-*.Rmd`.

The skill's directory name is the one name for the topic. Everything else is derived from it by `guide_index()` in `dev_scripts/build_guides.R` and must not be restated anywhere: the article is `vignettes/guide-<skill>.Rmd`, its title is `"Guide: "` followed by the skill's own H1, and other documents link to it as "the guide to \<H1 lowercased\>". So renaming a skill directory renames the article, and rewording its H1 rewords the title and every cross-reference to it. `dev_scripts/build_wiki.R` sources the same function rather than keeping its own copy of the mapping.

The one exception is `upgrade-mizer-code`, which is not a guide: it overrides all three names in the `guide_topics` table so that it stays `vignettes/upgrading.Rmd`, titled "Upgrading your mizer code". Nothing else should use those overrides. An article that renames because its skill is named differently is a rename like any other: add a `redirects:` entry in `pkgdown/_pkgdown.yml` and a row to the table in `inst/skills/upgrade-mizer-code/references/mizer-<version>.md`, rather than pinning the old name.

### Adding or modifying a skill

1. Edit the markdown file at `inst/skills/<topic>/SKILL.md`.
2. Wrap content meant strictly for AI agents (diagnostic decision trees, symptom tables) in `<!-- agent-only --> ... <!-- /agent-only -->`; the guide builder drops those blocks. Wrap content meant strictly for the article in `<!-- article-only --> ... <!-- /article-only -->`; the guide builder keeps it and drops the markers, while `mizerAgents::setup_mizer_agent()` drops the whole block as it installs the skill. Each side of the fence takes the half it wants and a topic still lives in one file.

   What belongs in an article-only block is a **demonstration**, not a definition: the chunk whose value is the output it produces — a plot, a printed range — which an agent cannot see and would only read past. The code being demonstrated stays in the skill body, where an agent asked to write something similar will actually find it; splitting the other way hides the function signatures an extension has to get right.

   Fence style is independent of this and purely mechanical, wherever either appears: ```` ```r ```` becomes `eval=FALSE` in the article, and ```` ```{r label} ```` passes through and *is* evaluated. So a definition the article's demonstrations depend on is written as an evaluated chunk in the skill body — that is what keeps the demonstration runnable without duplicating the definition into the block. Keep evaluated chunks cheap, and remember they have to keep working.

   The article-only stripping lives in mizerAgents (`.strip_article_only()`), so a mizer newer than the installed mizerAgents will hand an agent the blocks intact — noise, not breakage, and it clears when mizerAgents is updated.
3. Regenerate all guide vignettes:
   ```r
   source("dev_scripts/build_guides.R")
   build_guides()
   ```
4. If modifying behaviour or deprecating code, add the prose to `inst/skills/upgrade-mizer-code/references/mizer-<version>.md` and a symptom row to that skill's `SKILL.md` (together they regenerate `vignettes/upgrading.Rmd`).

---

## 2. The AI documentation index (`llms.txt`)

`llms.txt` follows the [llms.txt standard](https://llmstxt.org/): a concise package summary and categorized function index. It lives in two places — `docs/llms.txt`, served at <https://sizespectrum.org/mizer/llms.txt>, and `inst/llms.txt`, installed with the package so that `mizerAgents::setup_mizer_agent()` reads the index matching the installed mizer version.

`pkgdown::build_site()` writes a raw `docs/llms.txt` prefixed with `README.md` (badges and HTML). `build_llms()` replaces that preamble with `pkgdown/llms-header.md` — the file to edit to change the overview text — and syncs both locations:

```r
source("dev_scripts/build_llms.R")
build_llms()
```

The site build (section 3) already runs this. Call it by hand only after a direct `pkgdown::build_site()`, or to pick up an edit to `pkgdown/llms-header.md` without rebuilding the site.

`llms-full.txt` has been removed: function signatures are fetched dynamically from the installed package / R help system.

---

## 3. Building the website

Build the full site with the wrapper, not with `pkgdown::build_site()` on its own:

```bash
Rscript dev_scripts/build_site.R     # optional argument: the package root
```
```r
source("dev_scripts/build_site.R")
build_mizer_site()
```

It deletes any `vignettes/*.html`, runs `pkgdown::build_site()`, then runs `build_llms()` (section 2).

That deletion is the reason the wrapper exists. Knitting a vignette by hand — the RStudio Knit button, `rmarkdown::render()` — leaves a standalone HTML file beside its `.Rmd`, and pkgdown copies loose files from `vignettes/` into `docs/articles/`, so it lands on top of the formatted page for that article: the page silently loses its navbar, sidebar and styling, and nothing shows the regression until someone opens it on the website. The files are Git-ignored artefacts, so deleting them costs nothing. Any pkgdown call you make yourself skips this cleanup — delete stray HTML first.

### `_pkgdown.yml`

When adding, renaming or removing exported functions or vignettes, add them under `reference:` or `articles:` respectively.

### Fast iteration on specific pages

A full build is slow; rebuild single targets during development:

```r
pkgdown::build_article("model_description")
pkgdown::build_reference(topic = "getBiomass")
pkgdown::build_home()    # from index.md / README.md
pkgdown::build_news()
```

---

## 4. Synchronising navbar to `mizerCourse`

If `navbar:` in `_pkgdown.yml` changes, sync it with the mizer course website:

```bash
Rscript dev_scripts/sync_course_navbar.R [--course ../mizerCourse]
```

---

## 5. Developer wiki documentation

Developer skills in `.claude/skills/<topic>.md` are the single source of truth for AI agents and the GitHub wiki developer guides in `../mizer.wiki/`. Edit the skill, wrapping agent-only content as in section 1, then regenerate:

```bash
Rscript dev_scripts/build_wiki.R     # or source() it and call build_wiki()
```

This autolinks function names to their pkgdown reference pages and rewrites skill cross-references into wiki links.

---

## 6. Verification checklist before submitting a PR

1. `devtools::document()` runs cleanly with no Rd warnings.
2. `build_guides()` run if any `inst/skills/` file changed.
3. `build_wiki()` run if any `.claude/skills/` file changed.
4. `Rscript dev_scripts/build_site.R` run if any vignette, skill or `_pkgdown.yml` entry changed (it refreshes `llms.txt` too), and no `vignettes/*.html` left behind.
5. `_pkgdown.yml` updated for new exports or vignettes.
6. `Collate:` updated in `DESCRIPTION` if any new file was added in `R/`.
