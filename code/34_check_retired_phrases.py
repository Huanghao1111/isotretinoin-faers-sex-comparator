#!/usr/bin/env python
# ==============================================================================
# Retired phrases.
#
# Every number below was, at some point today, true in one file and false in
# another, or unsupported by any deposited table. Each was removed by hand. A
# hand removal lasts until the next edit; this list does not.
#
# When a phrase is retired, add it here with the reason. The check fails if it
# reappears in any document that could be uploaded.
# ==============================================================================
import re
import sys
import zipfile
from pathlib import Path

ROOT = Path("C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45")

RETIRED = [
    # (phrase, reason)
    ("1.26 where age is unrecorded", "no table holds this value; the two known-age strata are 1.28 and 1.51"),
    ("98% of isotretinoin", "no deposited table gives a lawyer share by event"),
    ("53 female versus 25 male", "TableS10 has no psychiatric-indication group"),
    ("TableS41", "the table is sheet S2 of the supplemental workbook now"),
    ("01b through", "the code runs to 32, not 20"),
    ("would be aimed at the wrong sex", "replaced: the key message names only a direction"),
    ("157 passed", "the verification suite prints 183"),
    ("29,844", "stale comparator removal count; the correct figure is 30,311"),
    ("73.5%", "the indication restriction removes 66.9%, not 73.5%"),
    ("9 to 15", "the sexual-term comparator factors run 4.7 to 14.5"),
    ("one thirty-fifth", "finasteride carries under 3% of the reports"),
    ("0.58 to 0.90", "the range across summaries and eras is 0.37 to 0.90"),
    ("a step added after the first complete analysis", "the disclosure now sits in the limitation paragraph"),
    ("never formally tested", "the cited meta-regression tested sex"),
    ("1.6 million users", "1,625,891 are participants, not users"),
    ("the reporting signal points to men", "the key message is now the property of the comparison"),
    ("remained after deduplication", "1,148,083 is what survives the indication restriction"),
]

# 99b_* files are dated reviews: they quote the text as it stood then, so they
# are allowed to contain retired wording. Everything else must be current.
FILES = [
    *(f for f in sorted((ROOT / "manuscript").glob("*.md")) if not f.name.startswith("99b_")),
    *sorted((ROOT / "supplement").glob("NOTE_*.md")),
    *sorted((ROOT / "supplement").glob("README.md")),
    ROOT / "README.md",
]

# the .docx files are what get uploaded, so they are scanned too
DOCX = sorted((ROOT / "submission").glob("*.docx"))


def docx_text(path):
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8")
    xml = re.sub(r"</w:p>", "\n", xml)
    return re.sub(r"<[^>]+>", "", xml)


def main():
    corpus = {}
    for f in FILES:
        if f.exists():
            corpus[f.name] = f.read_text(encoding="utf-8")
    for f in DOCX:
        corpus[f.name] = docx_text(f)

    hits = 0
    for phrase, reason in RETIRED:
        where = [n for n, t in corpus.items() if phrase in t]
        if where:
            hits += 1
            print(f'[FAIL] "{phrase}"')
            print(f"       still in: {', '.join(where)}")
            print(f"       why retired: {reason}")
        else:
            print(f'[PASS] "{phrase[:52]}" absent')

    print(f"\n===== {len(RETIRED) - hits} passed, {hits} failed =====")
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
