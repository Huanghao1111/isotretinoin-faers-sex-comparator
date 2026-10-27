#!/usr/bin/env python
# ==============================================================================
# Adds continuous line numbering to a .docx by injecting <w:lnNumType> into every
# section's properties. officer does not expose this and Word's own UI only sets
# it per section, so the package is edited directly. Idempotent: running it twice
# does not duplicate the element.
#
# Usage: python code/23b_line_numbers.py manuscript/xxx.docx
# ==============================================================================
import re
import shutil
import sys
import zipfile
from pathlib import Path

path = Path(sys.argv[1] if len(sys.argv) > 1 else
            "manuscript/JAAD_BriefReport_isotretinoin_sex.docx")
if not path.exists():
    sys.exit(f"missing {path}")

with zipfile.ZipFile(path) as z:
    names = z.namelist()
    blobs = {n: z.read(n) for n in names}

xml = blobs["word/document.xml"].decode("utf-8")

if "w:lnNumType" in xml:
    print("line numbers already present")
else:
    # sectPr is the last child of w:body and also appears inside w:pPr for
    # section breaks; the numbering element has a fixed position in the schema,
    # after pgSz/pgMar/paperSrc and before pgBorders/cols, so it is inserted
    # directly after the last margin element present.
    marker = "<w:lnNumType w:countBy=\"1\" w:restart=\"continuous\"/>"

    def add_ln(m):
        # The section properties are written by officer with self-closing margin
        # elements, so a closing-tag anchor never matches. Insert after the last
        # of the elements that must precede lnNumType in the schema.
        block = m.group(0)
        for pattern in (r"<w:pgMar\b[^>]*/>", r"<w:pgSz\b[^>]*/>",
                        r"<w:sectPr\b[^>]*>"):
            hit = re.search(pattern, block)
            if hit:
                return block[:hit.end()] + marker + block[hit.end():]
        return block

    new_xml, n = re.subn(r"<w:sectPr\b.*?</w:sectPr>", add_ln, xml, flags=re.S)
    if n == 0:
        sys.exit("no w:sectPr found: cannot add line numbers")
    blobs["word/document.xml"] = new_xml.encode("utf-8")
    print(f"line numbers added to {n} section(s)")

    tmp = path.with_suffix(".tmp.docx")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as out:
        for n_ in names:
            out.writestr(n_, blobs[n_])
    shutil.move(tmp, path)

with zipfile.ZipFile(path) as z:
    check = z.read("word/document.xml").decode("utf-8")
print("lnNumType count:", check.count("w:lnNumType"))
