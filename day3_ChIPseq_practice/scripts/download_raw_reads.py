#!/usr/bin/env python3
"""Download raw ChIP-seq reads for the course practice."""

from __future__ import annotations

import argparse
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin

import requests
from tqdm import tqdm


SAMPLES = {
    "MoPh7_H3K9me3": {
        "url": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-01-31_BGI/C-G/MoPh7_H3K9me3_r1/",
    },
    "MoPh7_H3K27Ac": {
        "url": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-01-31_BGI/C-G/MoPh7_H3K27Ac_r1/",
    },
    "MoPh7_H2A119Ub": {
        "url": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-01-31_BGI/C-G/MoPh7_H2A199Ub_r1/",
    },
    "MoPh7_input": {
        "url": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-01-31_BGI/C-G/MoPh7_H3K9me3_r1_i/",
    },
    "MoPh7_histone_input": {
        "url": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-01-31_BGI/C-G/MoPh7_H3K9me3_r1_i/",
    },
    "MoPh7_CTCF": {
        "url": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-12-15_BGI/MoPh7/",
    },
    "MoPh7_CTCF_input": {
        "url": "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-12-15_BGI/MoPh7_input/",
    },
}

DEFAULT_SAMPLES = [
    "MoPh7_H3K9me3",
    "MoPh7_H3K27Ac",
    "MoPh7_H2A119Ub",
    "MoPh7_input",
    "MoPh7_CTCF",
    "MoPh7_CTCF_input",
]


class LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        for key, value in attrs:
            if key == "href" and value:
                self.links.append(value)


def list_fastq_urls(base_url: str) -> list[str]:
    response = requests.get(base_url, timeout=60, verify=False)
    response.raise_for_status()

    parser = LinkParser()
    parser.feed(response.text)

    urls = []
    for link in parser.links:
        if re.search(r"\.f(ast)?q\.gz$", link):
            urls.append(urljoin(base_url, link))
    return sorted(urls)


def pick_mates(urls: list[str]) -> tuple[str, str]:
    r1 = [u for u in urls if re.search(r"(_1|_R1|R1_).*\.f(ast)?q\.gz$", Path(u).name)]
    r2 = [u for u in urls if re.search(r"(_2|_R2|R2_).*\.f(ast)?q\.gz$", Path(u).name)]
    if len(r1) != 1 or len(r2) != 1:
        names = "\n".join(Path(u).name for u in urls)
        raise RuntimeError(f"Could not detect exactly one R1 and one R2:\n{names}")
    return r1[0], r2[0]


def download(url: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"Skip existing: {out_path}")
        return

    with requests.get(url, stream=True, timeout=60, verify=False) as response:
        response.raise_for_status()
        total = int(response.headers.get("content-length", 0))
        with out_path.open("wb") as handle:
            progress = tqdm(
                response.iter_content(chunk_size=1024 * 1024),
                total=total // (1024 * 1024) if total else None,
                unit="MiB",
                desc=out_path.name,
            )
            for chunk in progress:
                if chunk:
                    handle.write(chunk)


def download_sample(sample: str, out_dir: Path, dry_run: bool = False) -> None:
    base_url = SAMPLES[sample]["url"]
    print(f"\n== {sample}")
    urls = list_fastq_urls(base_url)
    r1, r2 = pick_mates(urls)

    if dry_run:
        print(f"R1: {r1}")
        print(f"R2: {r2}")
        return

    sample_dir = out_dir / sample
    download(r1, sample_dir / f"{sample}_R1.fastq.gz")
    download(r2, sample_dir / f"{sample}_R2.fastq.gz")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample",
        action="append",
        choices=sorted(SAMPLES),
        help="Sample to download. Can be used several times. Default: all samples.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/raw"),
        help="Output directory. Default: data/raw",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print detected FASTQ URLs without downloading files.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    samples = args.sample or DEFAULT_SAMPLES
    requests.packages.urllib3.disable_warnings()  # self-signed certificate on storage
    for sample in samples:
        download_sample(sample, args.out_dir, dry_run=args.dry_run)
    print("\nDone.")


if __name__ == "__main__":
    main()
