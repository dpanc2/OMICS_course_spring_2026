#!/usr/bin/env python3

from __future__ import annotations

import ssl
from pathlib import Path
from urllib.parse import quote, urljoin
from urllib.request import Request, urlopen


BASE_URL = "https://genedev.bionet.nsc.ru/ftp/by_User/DashaPanchenko/OMICS_course_spring_2026/day4/bismark_output/MoPh7/"
DATA_DIR = Path("data/bismark")

FILES = [
    "MoPh7_1_bismark_bt2_PE_report.txt",
    "MoPh7_1_bismark_bt2_pe.CpG_report.txt",
    "MoPh7_1_bismark_bt2_pe.M-bias.txt",
    "MoPh7_1_bismark_bt2_pe.bedGraph.gz",
    "MoPh7_1_bismark_bt2_pe.bismark.cov.gz",
    "MoPh7_1_bismark_bt2_pe.cytosine_context_summary.txt",
    "MoPh7_1_bismark_bt2_pe_splitting_report.txt",
]


def download_file(filename: str) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    output_path = DATA_DIR / filename

    if output_path.exists() and output_path.stat().st_size > 0:
        print(f"Skip existing: {output_path}")
        return

    url = urljoin(BASE_URL, quote(filename))
    request = Request(url, headers={"User-Agent": "OMICS-course-downloader"})
    context = ssl._create_unverified_context()

    print(f"Download: {url}")
    with urlopen(request, context=context, timeout=120) as response:
        total = response.headers.get("Content-Length")
        total_text = f" ({int(total) / 1024 / 1024:.1f} MiB)" if total else ""
        print(f"Saving: {output_path}{total_text}")
        with output_path.open("wb") as out:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)


def main() -> None:
    for filename in FILES:
        download_file(filename)
    print("Done.")


if __name__ == "__main__":
    main()
