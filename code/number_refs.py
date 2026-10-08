"""Turn {key} placeholders into [n] and build the reference list.

Numbers are assigned by order of first appearance, which is the house rule and
also what the comparison manuscript does. Orphans (in the pool, never cited) and
unknown keys (cited, not in the pool) are both reported rather than silently
dropped, because a reference list that quietly loses an entry is worse than one
that fails loudly.

  python3 number_refs.py manuscript_src.md manuscript.md
  python3 number_refs.py manuscript_src.md manuscript.md --author-year

The second form emits (Author, Year) in the text and an alphabetical reference
list, for journals that require it. The house default is numbered by order of
appearance; only the output changes, never the source file.
"""
import re
import sys
from pathlib import Path

import refs

# a placeholder is {key} or {key1,key2}; the look-behind keeps subscripts such
# as p_{i} and superscripts such as e^{-kT} from being read as citations
PAT = re.compile(r"(?<![_^{])\{([a-zA-Z][a-zA-Z0-9_,\s]*)\}(?!\})")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    author_year = "--author-year" in sys.argv
    src = Path(args[0] if args else "manuscript_src.md")
    out = Path(args[1] if len(args) > 1 else "manuscript.md")
    text = src.read_text()

    order, unknown = [], []

    def repl(m):
        keys = [k.strip() for k in m.group(1).split(",")]
        nums = []
        for k in keys:
            if k not in refs.POOL:
                unknown.append(k)
                continue
            if k not in order:
                order.append(k)
            nums.append(order.index(k) + 1)
        if not nums:
            return m.group(0)
        # collapse runs of three or more into a range, as the comparison
        # manuscript does with [14-16]
        nums = sorted(set(nums))
        parts, i = [], 0
        while i < len(nums):
            j = i
            while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
                j += 1
            parts.append(str(nums[i]) if j - i < 2
                         else f"{nums[i]}-{nums[j]}")
            if j - i == 1:
                parts[-1] = f"{nums[i]},{nums[j]}"
            i = j + 1
        return "[" + ",".join(parts) + "]"

    def cite_author_year(m):
        keys = [k.strip() for k in m.group(1).split(",") if k.strip() in refs.POOL]
        for k in keys:
            if k not in order:
                order.append(k)
        parts = []
        for k in keys:
            e = refs.POOL[k]
            who = e.split(" (")[0]
            first = who.split(",")[0]
            multi = ("&" in who) or ("et al" in who)
            yr = re.search(r"\((\d{4}|n\.d\.)\)", e)
            parts.append(f"{first} et al., {yr.group(1)}" if multi
                         else f"{first}, {yr.group(1)}")
        return "(" + "; ".join(parts) + ")" if parts else m.group(0)

    body = PAT.sub(cite_author_year if author_year else repl, text)
    if author_year:
        reflist = "\n\n".join(refs.POOL[k] for k in
                               sorted(order, key=lambda k: refs.POOL[k].lower()))
    else:
        reflist = "\n\n".join(f"{i+1}. {refs.POOL[k]}"
                               for i, k in enumerate(order))
    body = body.replace("<<REFERENCES>>", reflist)

    # word count of the main text, Introduction to Discussion, for the header
    m = re.search(r"\nIntroduction\n(.*?)\n(Declaration|References)\n", body,
                  re.S)
    if m:
        words = len(re.sub(r"\[[\d,\-]+\]", "", m.group(1)).split())
        body = re.sub(r"(Word count \(main text, Introduction to Discussion\): )"
                      r"[\d,]+", lambda mm: mm.group(1) + f"{words:,}", body)
        body = re.sub(r"(Number of references: )\d+",
                      lambda mm: mm.group(1) + str(len(order)), body)

    # the abstract length on the title page is measured too. It was a typed
    # literal, and a typed literal goes stale the moment the abstract is cut
    # to a journal's limit.
    ma = re.search(r"\nAbstract\n(.*?)\nKeywords:", body, re.S)
    if ma:
        body = re.sub(r"(Abstract: )\d+( words)",
                      lambda mm: mm.group(1) + str(len(ma.group(1).split()))
                      + mm.group(2), body)
    out.write_text(body)

    # ---- report
    print(f"wrote {out}")
    print(f"references cited: {len(order)}; pool: {len(refs.POOL)}")
    orphans = [k for k in refs.POOL if k not in order]
    if orphans:
        print("ORPHANS (in pool, never cited):", ", ".join(orphans))
    if unknown:
        print("UNKNOWN KEYS (cited, not in pool):", ", ".join(sorted(set(unknown))))

    secs = ["Introduction", "Methods", "Results", "Discussion"]
    bounds = {s: body.find("\n" + s + "\n") for s in secs}
    # the end of the body is whichever of these comes first and exists: taking
    # the later one silently counted the whole Declaration block as Discussion
    ends = [body.find("\nDeclaration\n"), body.find("\nReferences\n")]
    bounds["end"] = min(x for x in ends if x > 0)
    keys = secs + ["end"]
    print(f"\n{'section':14s}{'words':>7}{'markers':>9}{'events':>8}"
          f"{'ev/1000w':>10}")
    tot_w = tot_e = 0
    for a, b in zip(keys, keys[1:]):
        seg = body[bounds[a]:bounds[b]]
        w = len(re.sub(r"\[[\d,\-]+\]", "", seg).split())
        mk = re.findall(r"\[[\d,\-]+\]", seg)
        ev = 0
        for x in mk:
            for part in x.strip("[]").split(","):
                ev += (int(part.split("-")[1]) - int(part.split("-")[0]) + 1
                       if "-" in part else 1)
        tot_w += w; tot_e += ev
        print(f"{a:14s}{w:7d}{len(mk):9d}{ev:8d}{1000*ev/max(w,1):10.1f}")
    print(f"{'BODY':14s}{tot_w:7d}{'':9s}{tot_e:8d}{1000*tot_e/max(tot_w,1):10.1f}")
    sents = re.split(r"(?<=[.!?])\s+", body[bounds["Introduction"]:bounds["end"]])
    wc = sum(1 for s in sents if re.search(r"\[\d", s))
    print(f"sentences {len(sents)}, with a citation {wc} "
          f"({100*wc/max(len(sents),1):.0f}%)")


if __name__ == "__main__":
    main()
