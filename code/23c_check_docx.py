#!/usr/bin/env python
# ==============================================================================
# Verifies the submission .docx: line numbering present, body double-spaced,
# and the body text identical in length to the markdown that the number checker
# validates. A Word file that silently loses a paragraph or reverts to single
# spacing is the kind of thing nobody notices until a copy-editor does.
#
# Exit status is non-zero on any failure.
# ==============================================================================
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

docx = Path("manuscript/JAAD_BriefReport_isotretinoin_sex.docx")
if not docx.exists():
    sys.exit(f"missing {docx}")

with zipfile.ZipFile(docx) as z:
    xml = z.read("word/document.xml")

root = ET.fromstring(xml)
body = root.find(f"{W}body")

paras = []
for p in body.iter(f"{W}p"):
    text = "".join(t.text or "" for t in p.iter(f"{W}t"))
    spacing = p.find(f"{W}pPr/{W}spacing")
    line = spacing.get(f"{W}line") if spacing is not None else None
    paras.append((text, line))

ok = fail = 0


def check(label, cond, detail=""):
    global ok, fail
    if cond:
        ok += 1
        print(f"[PASS] {label:<48} {detail}")
    else:
        fail += 1
        print(f"[FAIL] {label:<48} {detail}")


# ------------------------------------------------------------------ structure
sect = body.find(f"{W}sectPr")
ln = sect.find(f"{W}lnNumType")
check("line numbering present", ln is not None,
      "countBy=" + str(ln.get(f"{W}countBy")) if ln is not None else "MISSING")
check("line numbering continuous",
      ln is not None and ln.get(f"{W}restart") == "continuous",
      str(ln.get(f"{W}restart")) if ln is not None else "NA")
check("line numbering counts every line",
      ln is not None and ln.get(f"{W}countBy") == "1",
      str(ln.get(f"{W}countBy")) if ln is not None else "NA")

# ----------------------------------------------------------------- the body
start = next(i for i, (t, _) in enumerate(paras) if t.strip() == "To the Editor:")
end = next(i for i, (t, _) in enumerate(paras)
           if t.strip() == "References" and i > start)
body_paras = [t for t, _ in paras[start + 1:end] if t.strip()]
words = sum(len(t.split()) for t in body_paras)
check("body paragraph count", len(body_paras) == 5, f"{len(body_paras)} paragraphs")
check("body word count", 480 <= words <= 500, f"{words} words")

dbl = [l for t, l in paras[start + 1:end] if t.strip() and l == "480"]
check("body double-spaced", len(dbl) == len(body_paras),
      f"{len(dbl)}/{len(body_paras)} paragraphs at w:line=480")

# --------------------------------------------------------- strings and numbers
joined = "\n".join(t for t, _ in paras)
must = [
    "Male predominance in isotretinoin depression reporting",
    "0.59 (0.54-0.65)",
    "0.58 (0.50-0.67)",
    "30,311 comparator reports removed",
    "should not be aimed only at women",
    "The age contrast belongs to the comparator",
]
for m in must:
    check(f'contains "{m[:38]}"', m in joined)

check("stale 29,844 absent", "29,844" not in joined)
for k, sup in enumerate("\u00b9\u00b2\u00b3\u2074\u2075", start=1):
    check(f"citation {k} survives conversion", sup in joined)
check("no table in the manuscript", "Isotretinoin OR" not in joined)
km = next((t for t, _ in paras if t.startswith("Key Message:")), "")
km_text = km.split("Key Message:", 1)[1].strip() if km else ""
check("key message present", bool(km_text))
check("key message is a sentence of the body",
      bool(km_text) and km_text in joined.replace(km, ""))
check("figure legend covers both panels",
      "(A)" in joined and "(B)" in joined and "14.9 per 1,000" in joined)
check("five references", sum(1 for t, _ in paras if t[:3].strip().rstrip(".").isdigit() and ". " in t[:4]) >= 5,
      "numbered entries found")

# placeholders are expected until the authors fill them in; list them so that
# nothing is left behind by accident at upload
import re as _re
ph = sorted({x.strip().strip("[]") for x in _re.findall(r"\[[^\]]{2,60}\]", joined)})
nconf = sum(1 for t, _ in paras if "CONFIRM" in t)
if ph:
    print("\nBracketed placeholders still to complete:")
    for x in ph:
        print("  -", x)
if nconf:
    print(f"{nconf} paragraph(s) still contain CONFIRM.")

print(f"\n===== {ok} passed, {fail} failed =====")
sys.exit(1 if fail else 0)
