#!/usr/bin/env python
"""Build a self-contained HTML rendering of the astrovodka documentation.

The Markdown under ``docs/<lang>/`` is the single source of truth; this script
renders it to ``docs/html/<lang>/`` so the same documentation can be browsed
offline without a network connection or a Markdown viewer. The HTML is
committed, and ``tests/test_docs_sync.py`` fails if it drifts from the
Markdown, so the two cannot disagree the way the previously hand-maintained
HTML and Markdown trees did.

Deliberately uses only the standard library: no Markdown package, no build
tooling, no network. It implements the subset of Markdown these pages actually
use — ATX headings, pipe tables, fenced and indented code, bullet and ordered
lists, blockquotes, thematic breaks, and inline code/bold/italic/links.

Usage
-----
    python docs/build_html.py            # regenerate docs/html/
    python docs/build_html.py --check    # exit non-zero if it is out of date
"""

from __future__ import annotations

import argparse
import re
import sys
import unicodedata
from html import escape
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DOCS = REPO_ROOT / "docs"
OUT = DOCS / "html"
LANGS = ("en", "zh")

SITE_TITLE = "astrovodka documentation"

# Which language the reader is offered when they land on docs/html/index.html.
DEFAULT_LANG = "en"
LANG_LABEL = {"en": "English", "zh": "简体中文"}

# Links to files that live outside docs/ cannot be resolved from the rendered
# tree. Point them at the repository instead of emitting a dead relative link.
GITHUB_BLOB = "https://github.com/LRayleighJ/astrovodka/blob/master/"


# --------------------------------------------------------------------------
# Inline rendering
# --------------------------------------------------------------------------


def slugify(text: str) -> str:
    """Reproduce GitHub's heading-anchor algorithm.

    Lowercase, drop anything that is not alphanumeric/space/hyphen/underscore,
    then replace spaces with hyphens. Unicode letters are kept, so Chinese
    headings get usable anchors.
    """
    text = re.sub(r"<[^>]+>", "", text)
    text = text.strip().lower()
    text = "".join(
        ch
        for ch in text
        if ch == " " or ch == "-" or ch == "_" or unicodedata.category(ch)[0] in "LN"
    )
    return text.replace(" ", "-")


def inline(text: str, *, link_map=None) -> str:
    """Render inline Markdown to HTML.

    Code spans are extracted first so that Markdown syntax inside them is left
    alone, then links, bold and italic, and finally the code spans go back.
    """
    spans: list[str] = []

    def stash(match: re.Match) -> str:
        spans.append(match.group(1))
        return f"\x00{len(spans) - 1}\x00"

    text = re.sub(r"`([^`]+)`", stash, text)
    text = escape(text, quote=False)

    # Links, optionally with a title.
    def link(match: re.Match) -> str:
        label, href = match.group(1), match.group(2)
        if link_map is not None:
            href = link_map(href)
        return f'<a href="{escape(href, quote=True)}">{label}</a>'

    text = re.sub(r"\[([^\]]+)\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)", link, text)

    # Images are not used by these pages; render alt text rather than break.
    text = re.sub(r"!\[([^\]]*)\]\([^)]*\)", r"\1", text)

    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<![\w*])\*([^*\n]+)\*(?![\w*])", r"<em>\1</em>", text)

    def unstash(match: re.Match) -> str:
        return f"<code>{escape(spans[int(match.group(1))], quote=False)}</code>"

    return re.sub(r"\x00(\d+)\x00", unstash, text)


# --------------------------------------------------------------------------
# Block rendering
# --------------------------------------------------------------------------

_TABLE_SEP = re.compile(r"^\|?[\s:|-]+\|[\s:|-]*$")


def _split_row(row: str) -> list[str]:
    row = row.strip()
    if row.startswith("|"):
        row = row[1:]
    if row.endswith("|"):
        row = row[:-1]
    return [cell.strip() for cell in row.split("|")]


def render_table(lines: list[str], link_map) -> str:
    header = _split_row(lines[0])
    body = [_split_row(line) for line in lines[2:]]
    out = ['<table>', "<thead>", "<tr>"]
    out += [f"<th>{inline(c, link_map=link_map)}</th>" for c in header]
    out += ["</tr>", "</thead>", "<tbody>"]
    for row in body:
        out.append("<tr>")
        out += [f"<td>{inline(c, link_map=link_map)}</td>" for c in row]
        out.append("</tr>")
    out += ["</tbody>", "</table>"]
    return "\n".join(out)


_ITEM = re.compile(r"^(\s*)([-*+]|\d+[.)])\s+(.*)$")


def _render_list(block: list[str], link_map) -> str:
    """Render a (possibly nested) list block.

    Items are collected as ``(indent, ordered?, text)`` and folded into a tree
    by indentation before rendering, so a nested bullet list — the table of
    contents, for instance — produces properly nested ``<ul>`` elements inside
    its parent ``<li>`` instead of one flattened item.
    """
    items: list[tuple[int, bool, list[str]]] = []
    for raw in block:
        match = _ITEM.match(raw)
        if match:
            indent = len(match.group(1).expandtabs(4))
            ordered = match.group(2)[0].isdigit()
            items.append((indent, ordered, [match.group(3)]))
        elif items:
            # A wrapped continuation of the previous item.
            items[-1][2].append(raw.strip())
    if not items:
        return ""

    # ---- fold into a tree of nodes: {ordered, text, children} ----
    root: list[dict] = []
    # Stack of (indent, node-list-to-append-to). A deeper indent becomes a
    # child list of the most recent shallower item.
    stack: list[tuple[int, list]] = [(-1, root)]
    for indent, ordered, parts in items:
        while stack and indent <= stack[-1][0]:
            stack.pop()
        node = {"ordered": ordered, "text": " ".join(parts), "children": []}
        stack[-1][1].append(node)
        stack.append((indent, node["children"]))

    def render(nodes: list[dict]) -> list[str]:
        out: list[str] = []
        # A list is <ol> only if its first item is ordered, matching how
        # Markdown treats a change of marker as starting a new list.
        tag = "ol" if nodes and nodes[0]["ordered"] else "ul"
        out.append(f"<{tag}>")
        for node in nodes:
            body = inline(node["text"], link_map=link_map)
            if node["children"]:
                out.append(f"<li>{body}")
                out.extend(render(node["children"]))
                out.append("</li>")
            else:
                out.append(f"<li>{body}</li>")
        out.append(f"</{tag}>")
        return out

    return "\n".join(render(root))


def _is_table_start(lines: list[str], i: int) -> bool:
    return (
        i + 1 < len(lines)
        and "|" in lines[i]
        and _TABLE_SEP.match(lines[i + 1]) is not None
    )


def render_blocks(lines: list[str], link_map, anchors: list[tuple[int, str]]) -> str:
    """Render a block-level Markdown fragment to HTML."""
    out: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]

        # ---- fenced code ----
        fence = re.match(r"^\s*```(\w*)\s*$", line)
        if fence:
            lang = fence.group(1)
            i += 1
            buf = []
            while i < len(lines) and not re.match(r"^\s*```\s*$", lines[i]):
                buf.append(lines[i])
                i += 1
            i += 1
            cls = f' class="language-{lang}"' if lang else ""
            out.append(f"<pre><code{cls}>{escape(chr(10).join(buf))}</code></pre>")
            continue

        # ---- headings ----
        heading = re.match(r"^(#{1,6})\s+(.*)$", line)
        if heading:
            level = len(heading.group(1))
            text = heading.group(2).strip()
            anchor = slugify(text)
            anchors.append((level, anchor))
            out.append(
                f'<h{level} id="{escape(anchor, quote=True)}">'
                f"{inline(text, link_map=link_map)}</h{level}>"
            )
            i += 1
            continue

        # ---- thematic break ----
        if re.match(r"^\s*(-{3,}|\*{3,}|_{3,})\s*$", line):
            out.append("<hr>")
            i += 1
            continue

        # ---- table ----
        if _is_table_start(lines, i):
            block = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                block.append(lines[i])
                i += 1
            out.append(render_table(block, link_map))
            continue

        # ---- blockquote (may itself contain lists and code) ----
        if line.lstrip().startswith(">"):
            buf = []
            while i < len(lines) and (
                lines[i].lstrip().startswith(">") or lines[i].strip() == ""
            ):
                if lines[i].strip() == "" and not (
                    i + 1 < len(lines) and lines[i + 1].lstrip().startswith(">")
                ):
                    break
                buf.append(re.sub(r"^\s*>\s?", "", lines[i]))
                i += 1
            inner = render_blocks(buf, link_map, [])
            out.append(f"<blockquote>\n{inner}\n</blockquote>")
            continue

        # ---- lists ----
        if _ITEM.match(line):
            block = []
            # Consume the whole list, including nested items and wrapped
            # continuation lines, stopping at a blank line that is followed by
            # something which cannot be part of this list.
            while i < len(lines):
                cur = lines[i]
                if _ITEM.match(cur):
                    block.append(cur)
                    i += 1
                elif (
                    cur.strip()
                    and block
                    and not re.match(r"^\s*(#{1,6}\s|```|\||>)", cur)
                ):
                    block.append(cur)
                    i += 1
                elif (
                    cur.strip() == ""
                    and i + 1 < len(lines)
                    and (_ITEM.match(lines[i + 1]) or lines[i + 1].startswith("  "))
                ):
                    i += 1
                else:
                    break
            out.append(_render_list(block, link_map))
            continue

        # ---- blank ----
        if line.strip() == "":
            i += 1
            continue

        # ---- paragraph ----
        buf = [line]
        i += 1
        while (
            i < len(lines)
            and lines[i].strip() != ""
            and not re.match(r"^\s*(#{1,6}\s|```|>|\||[-*+]\s|\d+[.)]\s)", lines[i])
            and not re.match(r"^\s*(-{3,}|\*{3,}|_{3,})\s*$", lines[i])
        ):
            buf.append(lines[i])
            i += 1
        out.append(f"<p>{inline(' '.join(x.strip() for x in buf), link_map=link_map)}</p>")

    return "\n".join(out)


# --------------------------------------------------------------------------
# Page assembly
# --------------------------------------------------------------------------

CSS = """
:root{--fg:#1c1e21;--muted:#5c6570;--link:#0b62c4;--rule:#d8dee6;--bg:#fff;
--code-bg:#f4f6f8;--note-bg:#eef4fb;--head:#12395e}
@media (prefers-color-scheme:dark){:root{--fg:#e6e8eb;--muted:#9aa5b1;--link:#69a8f5;
--rule:#39404a;--bg:#15181c;--code-bg:#1e2227;--note-bg:#1b2733;--head:#9cc6f0}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);line-height:1.65;
font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",
"Noto Sans CJK SC","Noto Sans CJK JP","PingFang SC","Microsoft YaHei",sans-serif}
.wrap{max-width:900px;margin:0 auto;padding:2rem 1.25rem 4rem}
nav.lang{font-size:.85rem;color:var(--muted);margin-bottom:2rem;
display:flex;gap:1rem;flex-wrap:wrap;align-items:center}
nav.lang a{color:var(--link);text-decoration:none}
nav.lang a:hover{text-decoration:underline}
h1,h2,h3,h4,h5,h6{color:var(--head);line-height:1.3;margin:2rem 0 .75rem}
h1{font-size:1.9rem;border-bottom:2px solid var(--rule);padding-bottom:.4rem}
h2{font-size:1.45rem;border-bottom:1px solid var(--rule);padding-bottom:.3rem}
h3{font-size:1.15rem}h4{font-size:1rem}
a{color:var(--link)}
p,li{margin:.5rem 0}
code{background:var(--code-bg);padding:.12em .38em;border-radius:4px;
font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.9em}
pre{background:var(--code-bg);padding:1rem;border-radius:8px;overflow-x:auto;
border:1px solid var(--rule)}
pre code{background:none;padding:0;font-size:.86rem;line-height:1.5}
table{border-collapse:collapse;width:100%;margin:1rem 0;display:block;overflow-x:auto}
th,td{border:1px solid var(--rule);padding:.45rem .7rem;text-align:left;
vertical-align:top;font-size:.92rem}
th{background:var(--code-bg);font-weight:600}
blockquote{margin:1rem 0;padding:.6rem 1rem;border-left:4px solid var(--link);
background:var(--note-bg);border-radius:0 6px 6px 0}
blockquote p{margin:.3rem 0}
hr{border:0;border-top:1px solid var(--rule);margin:2.5rem 0}
footer{margin-top:4rem;padding-top:1rem;border-top:1px solid var(--rule);
color:var(--muted);font-size:.82rem}
.idx{display:grid;grid-template-columns:repeat(auto-fill,minmax(220px,1fr));gap:.6rem;
margin:1.5rem 0}
.idx a{display:block;padding:.7rem .9rem;border:1px solid var(--rule);border-radius:8px;
text-decoration:none;background:var(--code-bg)}
.idx a:hover{border-color:var(--link)}
"""


def page(title: str, body: str, lang: str, rel_root: str, extra_nav: str = "") -> str:
    """Wrap rendered body HTML in a complete document.

    ``rel_root`` is the relative path back to ``docs/html/``, so the language
    switch always lands on that language's landing page regardless of how deep
    the current page sits.
    """
    other = "zh" if lang == "en" else "en"
    switch = f'<a href="{rel_root}{other}/index.html">{LANG_LABEL[other]}</a>'
    return f"""<!DOCTYPE html>
<html lang="{lang}">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)}</title>
<style>{CSS}</style>
</head>
<body>
<div class="wrap">
<nav class="lang"><a href="{rel_root}index.html">{SITE_TITLE}</a>{switch}{extra_nav}</nav>
{body}
<footer>astrovodka &mdash; mock resolved HI observations of IllustrisTNG galaxies.
Generated from the Markdown sources by <code>docs/build_html.py</code>.</footer>
</div>
</body>
</html>
"""


def collect_anchors(lines: list[str]) -> list[tuple[int, str]]:
    anchors: list[tuple[int, str]] = []
    render_blocks(lines, None, anchors)
    return anchors


def make_link_map(src: Path, lang: str, known: set[Path]):
    """Return a function rewriting a Markdown link for the HTML tree."""

    def rewrite(href: str) -> str:
        if re.match(r"^(https?:|mailto:|#)", href):
            return href
        if href.startswith("$") or href.startswith("<"):
            return href
        target = (src.parent / href).resolve()
        try:
            rel = target.relative_to(DOCS)
        except ValueError:
            # Outside docs/ (a README, a source file, an example): send the
            # reader to the repository, where the relative path does resolve.
            return GITHUB_BLOB + href.lstrip("./")
        if href.endswith(".md"):
            if target not in known:
                return GITHUB_BLOB + href.lstrip("./")
            return str(rel.with_suffix(".html"))
        return str(rel)

    return rewrite


def build(check_only: bool = False) -> int:
    """Render every Markdown page. Returns the number of files written."""
    sources = {
        p.resolve(): p for p in sorted(DOCS.glob("*/*.md")) if p.parent.name in LANGS
    }
    sources.update(
        {p.resolve(): p for p in sorted(DOCS.glob("*/*/*.md")) if p.parent.parent.name in LANGS}
    )

    stale: list[str] = []
    written = 0

    for src in sources.values():
        lang = src.relative_to(DOCS).parts[0]
        rel = src.relative_to(DOCS / lang)
        target = OUT / lang / rel.with_suffix(".html")
        depth = len(rel.parts) - 1
        rel_root = "../" * (depth + 1)

        text = src.read_text(encoding="utf-8")
        lines = text.splitlines()
        body = render_blocks(lines, make_link_map(src, lang, set(sources)), [])
        title = next(
            (ln.lstrip("# ").strip() for ln in lines if ln.startswith("# ")),
            src.stem,
        )
        html = page(title, body, lang, rel_root)

        if check_only:
            if not target.is_file() or target.read_text(encoding="utf-8") != html:
                stale.append(str(target.relative_to(REPO_ROOT)))
            continue

        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(html, encoding="utf-8")
        written += 1

    # Language chooser at the top of the tree.
    if not check_only:
        entries = "\n".join(
            f'<a href="{lg}/index.html">{LANG_LABEL[lg]} &mdash; astrovodka</a>'
            for lg in LANGS
        )
        (OUT / "index.html").write_text(
            page(
                SITE_TITLE,
                f'<h1>{SITE_TITLE}</h1><div class="idx">{entries}</div>',
                DEFAULT_LANG,
                "",
            ),
            encoding="utf-8",
        )

    if check_only:
        if stale:
            print("docs/html/ is out of date. Regenerate with:")
            print("    python docs/build_html.py")
            for path in stale:
                print(f"  stale: {path}")
            return 1
        print(f"docs/html/ is up to date ({len(sources)} pages).")
        return 0

    print(f"Wrote {written} pages to {OUT.relative_to(REPO_ROOT)}/")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check",
        action="store_true",
        help="do not write; exit non-zero if docs/html/ is out of date",
    )
    return build(check_only=parser.parse_args().check)


if __name__ == "__main__":
    sys.exit(main())
