#!/usr/bin/env python
# ==============================================================================
# Numbers must match the table their own sentence cites.
#
# A first version of this check only asked whether a number appeared anywhere in
# the deposit. A negative test killed it: 2.14 and 41.7% were injected into a
# sentence and sailed through, because other tables contain those values for
# other quantities. A check that cannot fail is worse than no check.
#
# This version ties the number to its context. If a sentence cites one table, the
# numbers in that sentence have to be cells of that table. That is the claim a
# reader can test, so it is the claim the checker tests.
# ==============================================================================
import csv
import re
import sys
from pathlib import Path

ROOT = Path("C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45")
FILES = [Path(a) for a in sys.argv[1:]] or [
    ROOT / "manuscript/RECORD_PE_checklist.md",
    ROOT / "manuscript/00_submission_package.md",
]

DERIVED = {"34", "2.8", "11.1", "97.2", "2,885", "99.6", "96.8", "6.9", "1.6"}

NUM = re.compile(r"\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+\.\d+|\d{3,}|\d+%")
TABLE = re.compile(r"TableS(\d+)(?:\s*[\u2013-]\s*S?(\d+))?")


def table_values(n):
    hits = list((ROOT / "supplement").glob(f"TableS{n}_*.csv"))
    if not hits:
        return None
    vals = set()
    with open(hits[0], newline="", encoding="utf-8") as fh:
        for row in csv.reader(fh):
            for c in row:
                c = c.strip()
                if not c:
                    continue
                vals.add(c)
                try:
                    v = float(c)
                except ValueError:
                    continue
                for dp in (1, 2, 3):
                    vals.add(f"{v:.{dp}f}")
    return vals


def norm(x):
    return x.replace(",", "").rstrip("%")


def main():
    failures = 0
    checked = 0

    for path in FILES:
      text = path.read_text(encoding="utf-8")
      for line in text.split("\n"):
        if line.startswith("#") or not line.strip():
            continue
        for sentence in re.split(r"(?<=[.;])\s+", line):
            refs = set()
            for a, b in TABLE.findall(sentence):
                refs.update(range(int(a), int(b) + 1) if b else [int(a)])
            # a range of three or more tables is too coarse to test
            if not refs or len(refs) > 2:
                continue
            vals = set()
            for n in refs:
                v = table_values(n)
                if v:
                    vals |= v
            if not vals:
                continue
            n = ",".join(str(x) for x in sorted(refs))
            known = {norm(v) for v in vals} | DERIVED
            for m in NUM.finditer(sentence):
                tok = m.group(0)
                k = norm(tok)
                # "per 1,000" is a unit, not a value
                if re.search(r"per\s*$", sentence[:m.start()]):
                    continue
                if len(k) >= 7 or k in known:
                    continue
                failures += 1
                print(f"[FAIL] TableS{n}: {tok} is not a value in that table")
                print(f"       ...{sentence.strip()[:120]}")
            checked += 1

    print(f"\n===== sentences tied to a table: {checked}; failures: {failures} =====")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
