#!/usr/bin/env python3
"""
Assemble docs/llms-full.txt from docs/llms.txt plus inlined reference pages
and articles. Also copies llms-full.txt and llms.txt into ../mizerAGENTS/inst/.

Run from the mizer repo root:
    python3 dev_scripts/build_llms_full.py

Requires ../mizerAGENTS/inst/ to exist (the mizerAGENTS package repo must be
checked out as a sibling directory).

After running pkgdown::build_site() you should also restore the hand-edited
header in docs/llms.txt before running this script — see the wiki page
"Editing the Mizer Website" for instructions.
"""
import re
import os

DOCS = "docs"

# Reference pages to inline, in order of importance to a typical user.
REF_PAGES = [
    "MizerParams-class.md",
    "MizerSim-class.md",
    "species_params.md",
    "newMultispeciesParams.md",
    "project.md",
    "steady.md",
    "setFishing.md",
    "setResource.md",
    "setBevertonHolt.md",
    "getRates.md",
    "matchGrowth.md",
]

# Articles to inline (bibliography and image lines are stripped).
ARTICLES = [
    "model_description.md",
    "extending-mizer.md",
    "cheatsheet-analysis-and-plotting.md",
]


def read(path):
    with open(path) as f:
        return f.read()


def strip_images(text):
    return "\n".join(
        line for line in text.split("\n")
        if not re.match(r"^\s*!\[", line)
    )


def strip_references(text):
    m = re.search(r"\n## References\n", text)
    return text[: m.start()] if m else text


def process_ref(name):
    path = os.path.join(DOCS, "reference", name)
    return strip_images(read(path)).rstrip()


def process_article(name):
    path = os.path.join(DOCS, "articles", name)
    text = read(path)
    return strip_references(strip_images(text)).rstrip()


def main():
    parts = []

    parts.append(read(os.path.join(DOCS, "llms.txt")).rstrip())

    parts.append(
        "\n\n---\n\n"
        "# Full Reference Documentation\n\n"
        "The sections below contain the complete documentation for the most "
        "important functions and articles. They are included here so that an "
        "AI agent can answer questions about mizer without fetching individual "
        "pages.\n"
    )

    parts.append("\n## Reference pages\n")
    for name in REF_PAGES:
        path = os.path.join(DOCS, "reference", name)
        if os.path.exists(path):
            parts.append("\n---\n\n" + process_ref(name) + "\n")
        else:
            print(f"WARNING: not found — {path}")

    parts.append("\n## Key articles\n")
    for name in ARTICLES:
        path = os.path.join(DOCS, "articles", name)
        if os.path.exists(path):
            parts.append("\n---\n\n" + process_article(name) + "\n")
        else:
            print(f"WARNING: not found — {path}")

    output = "\n".join(parts)

    MIZER_AGENTS_INST = os.path.join("..", "mizerAGENTS", "inst")

    for outpath in [os.path.join(DOCS, "llms-full.txt"),
                    os.path.join(MIZER_AGENTS_INST, "llms-full.txt")]:
        with open(outpath, "w") as f:
            f.write(output)

    # Also bundle the lightweight index
    llms_index = read(os.path.join(DOCS, "llms.txt"))
    with open(os.path.join(MIZER_AGENTS_INST, "llms.txt"), "w") as f:
        f.write(llms_index)

    line_count = output.count("\n") + 1
    print(f"Written {line_count} lines ({len(output):,} bytes) to "
          f"docs/llms-full.txt, ../mizerAGENTS/inst/llms-full.txt, and "
          f"../mizerAGENTS/inst/llms.txt")


if __name__ == "__main__":
    main()
