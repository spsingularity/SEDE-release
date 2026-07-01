#!/usr/bin/env python3
"""
Assemble the build-time markdown: (1) convert informal inline math via mdmath, then
(2) lift the H1 title + the Abstract section into YAML metadata and add the author, so the
PDF renders a proper title block (title / author / date / abstract) and --toc lands AFTER
the abstract — the JCAP/PRD convention. SEDE.md itself is left untouched.
Run: python makedoc.py SEDE.md OUT.md
"""
import re, sys, datetime, mdmath

def main(inp, outp):
    raw = open(inp, encoding="utf-8").read()
    # 1. wrap inline math
    body, spans = mdmath.protect(raw)
    text = mdmath.restore(mdmath.RUN.sub(mdmath.convert_run, body), spans)
    lines = text.split("\n")

    # 2a. title = first '# ' heading
    title = ""
    ti = next(i for i, L in enumerate(lines) if L.startswith("# "))
    title = lines[ti][2:].strip()

    # 2b. abstract = body of the '## Abstract' section, up to the next '## ' or a '---' rule
    ai = next(i for i, L in enumerate(lines) if re.match(r"##\s+Abstract\s*$", L))
    j = ai + 1
    abs = []
    while j < len(lines) and not (lines[j].startswith("## ") or lines[j].strip() == "---"):
        abs.append(lines[j]); j += 1
    abstract = "\n".join(abs).strip()

    # 2c. body = everything from the first numbered section onward
    bi = next(i for i, L in enumerate(lines) if re.match(r"##\s+1\.", L))
    rest = "\n".join(lines[bi:])

    def yml(s):  # escape for a double-quoted YAML scalar
        return s.replace("\\", "\\\\").replace('"', '\\"')

    date = datetime.date.today().strftime("%-d %B %Y")
    abs_block = "\n".join(("  " + L) if L.strip() else "" for L in abstract.split("\n"))
    front = (
        "---\n"
        f'title: "{yml(title)}"\n'
        'author: "Stilian Pandev"\n'
        f'date: "{yml(date)}"\n'
        "abstract: |\n"
        f"{abs_block}\n"
        "---\n\n"
    )
    open(outp, "w", encoding="utf-8").write(front + rest + "\n")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
