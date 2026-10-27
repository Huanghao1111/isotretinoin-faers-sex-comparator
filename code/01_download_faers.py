"""
Download FAERS quarterly ASCII files from FDA.
Only quarters from 2012Q4 onward are available at fis.fda.gov/content/Exports/.
"""
import os
import sys
import time
import urllib.request
import traceback

RAW = r"C:\Users\lenovo\WorkBuddy\2026-10-17-00-49-45\data\raw"
os.makedirs(RAW, exist_ok=True)

BASE = "https://fis.fda.gov/content/Exports/faers_ascii_{y}q{q}.zip"
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"

MIN_SIZE = 500_000  # anything smaller is an error page


def quarters(start=(2012, 4), end=(2026, 4)):
    y, q = start
    while (y, q) <= end:
        yield y, q
        q += 1
        if q > 4:
            q = 1
            y += 1


def fetch(url, dest, tries=3):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=180) as r:
                data = r.read()
            if len(data) < MIN_SIZE:
                return None, f"too small ({len(data)}B)"
            with open(dest, "wb") as f:
                f.write(data)
            return len(data), "ok"
        except Exception as e:
            if attempt == tries - 1:
                return None, f"fail: {e}"
            time.sleep(5)
    return None, "fail"


def main():
    ok, miss, skip = 0, [], 0
    for y, q in quarters():
        name = f"faers_ascii_{y}q{q}.zip"
        dest = os.path.join(RAW, name)
        if os.path.exists(dest) and os.path.getsize(dest) >= MIN_SIZE:
            skip += 1
            continue
        size, status = fetch(BASE.format(y=y, q=q), dest)
        tag = f"{y}q{q}"
        if size:
            ok += 1
            print(f"[OK]   {tag}  {size/1e6:.1f} MB", flush=True)
        else:
            if os.path.exists(dest):
                os.remove(dest)
            if "too small" in str(status):
                miss.append(tag)
            print(f"[MISS] {tag}  {status}", flush=True)
    print("\n==== SUMMARY ====")
    print(f"downloaded={ok}  skipped(already)={skip}  unavailable={len(miss)}")
    print("unavailable:", " ".join(miss) if miss else "none")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
