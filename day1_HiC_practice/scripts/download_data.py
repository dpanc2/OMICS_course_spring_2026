#!/usr/bin/env python3
"""Download Hi-C maps, convert .hic to .mcool, and balance resolutions."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path
from urllib.parse import urlparse

import cooler
import requests
from tqdm import tqdm

try:
    from hic2cool import hic2cool_convert
except ImportError:
    from hic2cool.hic2cool_utils import hic2cool_convert


def filename_from_url(url: str) -> str:
    return Path(urlparse(url).path).name


def download_file(url: str, dest_dir: Path) -> Path:
    """Download a file from a URL unless it already exists."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    filepath = dest_dir / filename_from_url(url)

    if filepath.exists():
        print(f"{filepath.name} already exists. Skipping download.")
        return filepath

    print(f"Downloading {filepath.name}...")
    response = requests.get(url, stream=True, timeout=60)
    response.raise_for_status()

    total_size = int(response.headers.get("content-length", 0))
    block_size = 1024 * 1024

    with filepath.open("wb") as file_handle, tqdm(
        desc=filepath.name,
        total=total_size,
        unit="B",
        unit_scale=True,
        unit_divisor=1024,
    ) as progress:
        for chunk in response.iter_content(block_size):
            if chunk:
                progress.update(len(chunk))
                file_handle.write(chunk)

    return filepath


def prepare_local_hic(hic_path: Path, dest_dir: Path) -> Path:
    """Copy a local .hic into the data directory when needed."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    output_path = dest_dir / hic_path.name

    if hic_path.resolve() == output_path.resolve():
        return output_path

    if output_path.exists():
        print(f"{output_path.name} already exists. Skipping copy.")
        return output_path

    print(f"Copying {hic_path} -> {output_path}")
    shutil.copy2(hic_path, output_path)
    return output_path


def convert_hic_to_mcool(hic_path: Path, nproc: int = 4, force: bool = False) -> Path:
    """Convert a .hic file into a multi-resolution .mcool file."""
    mcool_path = hic_path.with_suffix(".mcool")

    if mcool_path.exists() and not force:
        print(f"{mcool_path.name} already exists. Skipping conversion.")
        return mcool_path

    if mcool_path.exists() and force:
        mcool_path.unlink()

    print(f"Converting {hic_path} -> {mcool_path}")
    hic2cool_convert(
        str(hic_path),
        str(mcool_path),
        resolution=0,
        nproc=nproc,
    )
    return mcool_path


def balance_mcool(mcool_path: Path, nproc: int = 4) -> None:
    """Balance all resolutions inside an .mcool file."""
    print(f"Balancing {mcool_path}...")
    resolutions = cooler.fileops.list_coolers(str(mcool_path))

    for resolution_group in resolutions:
        uri = f"{mcool_path}::{resolution_group}"
        print(f"  Balancing {uri}")
        cmd = [
            "cooler",
            "balance",
            "--force",
            "--ignore-diags",
            "2",
            "--mad-max",
            "5",
            "--min-nnz",
            "10",
            "--nproc",
            str(nproc),
            uri,
        ]
        subprocess.run(cmd, check=True)


def process_hic(hic_path: Path, nproc: int, force_convert: bool, skip_balance: bool) -> Path:
    mcool_path = convert_hic_to_mcool(hic_path, nproc=nproc, force=force_convert)
    if not skip_balance:
        balance_mcool(mcool_path, nproc=nproc)
    return mcool_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download Hi-C .hic maps, convert them to .mcool, and balance all resolutions."
    )
    parser.add_argument(
        "--url",
        nargs="+",
        default=None,
        help="One or more .hic URLs.",
    )
    parser.add_argument(
        "--hic",
        nargs="+",
        type=Path,
        default=[],
        help="One or more local .hic files to copy into data-dir and process.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data"),
        help="Directory for downloaded/copied .hic and generated .mcool files.",
    )
    parser.add_argument("--nproc", type=int, default=4, help="Number of processes for conversion/balancing.")
    parser.add_argument("--force-convert", action="store_true", help="Recreate .mcool files if they already exist.")
    parser.add_argument("--skip-balance", action="store_true", help="Convert only, without cooler balance.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    urls = args.url or []

    hic_files = []

    for url in urls:
        hic_files.append(download_file(url, args.data_dir))

    for hic_path in args.hic:
        hic_files.append(prepare_local_hic(hic_path, args.data_dir))

    if not hic_files:
        raise SystemExit("No input maps. Provide --url or --hic.")

    for hic_path in hic_files:
        process_hic(
            hic_path,
            nproc=args.nproc,
            force_convert=args.force_convert,
            skip_balance=args.skip_balance,
        )

    print("All Hi-C maps are ready.")


if __name__ == "__main__":
    main()
