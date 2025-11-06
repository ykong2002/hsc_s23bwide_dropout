#!/usr/bin/env python3
# fetch_all.py
#  1) Build the HSC coord-list for downloadCutout.py
#  2) (OPTIONAL) Run downloadCutout.py on that list (use --run-download)
#  3) Optionally fetch LoTSS cutouts
#  4) Optionally fetch VLASS cutouts from CADC

import os
import sys
import argparse
import subprocess
from collections import defaultdict
from urllib.parse import quote
import requests
import numpy as np
import pandas as pd

# -------------------------
# HSC list builder
# -------------------------

def build_hsc_cutout_list(
    survey_name: str,
    dropout: str,
    root: str,
    *,
    catalog_path: str = None,
    cutout_sizes=(2, 4, 8, 12, 30, 90, 150),  # arcsec (semi-size; sw/sh)
    filters=('g', 'r', 'i', 'z', 'y'),
    rerun='s23b_wide',
    output_txt: str = None,
) -> str:
    """
    Create a coord list text file for downloadCutout.py.

    The file is written to: {root}/{survey_name}_{dropout}.txt (by default)

    Columns:
      #? rerun filter ra dec sw sh name
    where `name` is a file path *without* ".fits" — downloadCutout.py will append.

    HSC files will save under:
      {root}/{survey_name}/{survey_name}_{dropout}_hsc/{tract}_{filter}_{size}.fits
    """
    # Resolve defaults (all under root)
    if catalog_path is None:
        catalog_path = os.path.join(root, f"{survey_name}/{survey_name}_{dropout}_dropout.csv")
    if output_txt is None:
        output_txt = os.path.join(root, f"{survey_name}_{dropout}.txt")

    # Load catalog
    df = pd.read_csv(catalog_path)

    # Flexible RA/Dec column names
    ra_col  = 'ra'  if 'ra'  in df.columns else 'RA'
    dec_col = 'dec' if 'dec' in df.columns else 'DEC'
    tract_col = 'tract'
    for col in (ra_col, dec_col, tract_col):
        if col not in df.columns:
            raise ValueError(f"Catalog must include a '{col}' column (have: {list[df.columns]})")

    header = "#? rerun    filter    ra        dec        sw      sh     name"
    out_lines = [header]

    # Ensure HSC output directory exists
    hsc_dir = os.path.join(root, f"{survey_name}", f"{survey_name}_{dropout}_hsc")
    os.makedirs(hsc_dir, exist_ok=True)

    # Deduplicate tract names if repeated
    tract_count = defaultdict(int)

    for _, row in df.iterrows():
        ra  = float(row[ra_col])
        dec = float(row[dec_col])
        tract_base = str(row[tract_col])

        n = tract_count[tract_base]
        tract_count[tract_base] += 1
        tract_tag = tract_base if n == 0 else f"{tract_base}_{n}"

        for cut in cutout_sizes:
            sw = f"{cut}asec"
            sh = f"{cut}asec"
            for f in filters:
                # Where downloadCutout.py will save: {name}.fits
                name_no_ext = os.path.join(hsc_dir, f"{tract_tag}_{f}_{cut}")
                line = f"{rerun:10} {f:7} {ra:.7f} {dec:.7f} {sw:6} {sh:6} {name_no_ext}"
                out_lines.append(line)

    # Write the list
    os.makedirs(root, exist_ok=True)
    with open(output_txt, "w") as fo:
        fo.write("\n".join(out_lines) + "\n")

    print(f"[HSC] Wrote list: {output_txt}  (lines: {len(out_lines)-1})")
    return output_txt


# -------------------------
# LoTSS cutout fetcher
# -------------------------

def fetch_lotss_cutouts(
    survey_name: str,
    dropout: str,
    catalog_path: str,
    radio_sizes_arcmin=(1, 3, 5),
    output_dir: str = None,
):
    df = pd.read_csv(catalog_path)
    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(catalog_path), f"../{survey_name}_{dropout}")
        output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # RA/Dec + tract columns
    ra_col  = 'RA'  if 'RA'  in df.columns else 'ra'
    dec_col = 'DEC' if 'DEC' in df.columns else 'dec'
    tract_col = 'tract'
    for col in (ra_col, dec_col, tract_col):
        if col not in df.columns:
            raise ValueError(f"LoTSS: catalog must include '{col}' (found {list(df.columns)})")

    tract_count = defaultdict(int)
    total = 0
    ok = 0

    for ra, dec, tract_base in zip(df[ra_col], df[dec_col], df[tract_col]):
        n = tract_count[tract_base]
        tract_count[tract_base] += 1
        tract_tag = tract_base if n == 0 else f"{tract_base}_{n}"
        pos = f"{ra} {dec}"
        encoded_pos = quote(pos)

        for size in radio_sizes_arcmin:
            total += 1
            url = f"https://lofar-surveys.org/dr2-cutout.fits?pos={encoded_pos}&size={size}"
            out_fn = os.path.join(output_dir, f"tract_{tract_tag}_lotss_{size}.fits")

            try:
                r = requests.get(url, timeout=120)
                if r.status_code == 200 and r.content:
                    with open(out_fn, "wb") as f:
                        f.write(r.content)
                    ok += 1
                    print(f"[LoTSS] saved: {out_fn}")
                else:
                    print(f"[LoTSS] failed ({r.status_code}) for {tract_tag} size={size}")
            except Exception as e:
                print(f"[LoTSS] error for {tract_tag} size={size}: {e}")

    print(f"[LoTSS] done: {ok}/{total} files saved to {output_dir}")


# -------------------------
# VLASS cutout fetcher (CADC)
# -------------------------

def fetch_vlass_cutouts_cadc(
    catalog_path: str,
    out_dir: str,
    *,
    radio_sizes_arcmin=(1, 3, 5),
    require_astroquery=True,
):
    """
    Minimal VLASS cutout fetcher using CADC:
      - query_region at the requested radius
      - get_image_list to obtain cutout URLs
      - download the last URL returned
    Saves to: {out_dir}/tract_{tract}_vlass_{size}.fits
    """
    # Imports
    if require_astroquery:
        try:
            from astroquery.cadc import Cadc
            from astropy import units as u
        except Exception as e:
            print(f"[VLASS] astroquery not available: {e} — skipping VLASS.")
            return

    # Load catalog
    df = pd.read_csv(catalog_path)

    os.makedirs(out_dir, exist_ok=True)
    cadc = Cadc()

    # Columns
    ra_col  = 'RA'  if 'RA'  in df.columns else 'ra'
    dec_col = 'DEC' if 'DEC' in df.columns else 'dec'
    tract_col = 'tract'
    for col in (ra_col, dec_col, tract_col):
        if col not in df.columns:
            raise ValueError(f"VLASS: catalog must include '{col}' (found {list(df.columns)})")

    # Deduplicate tract names
    tract_count = defaultdict(int)

    ok = 0
    total = 0

    for _, row in df.iterrows():
        ra = row[ra_col]
        dec = row[dec_col]
        tract_base = str(row[tract_col])

        # unique tract tag
        n = tract_count[tract_base]
        tract_count[tract_base] += 1
        tract = tract_base if n == 0 else f"{tract_base}_{n}"

        coords = f"{ra}, {dec}"

        for size in radio_sizes_arcmin:
            total += 1
            try:
                results = cadc.query_region(coords, radius=size * u.arcmin, collection='VLASS')
                image_list = cadc.get_image_list(results, coords, radius=size * u.arcmin)
            except Exception as e:
                print(f"[VLASS] CADC query failed for {tract} size={size}: {e}")
                continue

            if not image_list:
                print(f"[VLASS] no cutout URLs for {tract} (size={size}')")
                continue

            url = image_list[-1]  # simplest choice: last URL
            out_fn = os.path.join(out_dir, f"tract_{tract}_vlass_{size:g}.fits")

            print(f"[VLASS] Fetching {tract} {size} arcmin cutout…")
            try:
                r = requests.get(url, timeout=300)
                if r.status_code == 200 and r.content:
                    with open(out_fn, "wb") as f:
                        f.write(r.content)
                    ok += 1
                    print(f"[VLASS] Saved: {out_fn}")
                else:
                    print(f"[VLASS] Failed for {tract} size={size}: HTTP {r.status_code}")
            except Exception as e:
                print(f"[VLASS] Download error for {tract} size={size}: {e}")

    print(f"[VLASS] done: {ok}/{total} attempts wrote files to {out_dir}")



# -------------------------
# Utility: run downloadCutout.py
# -------------------------

def run_download_cutout(list_path: str, script_path: str = "downloadCutout.py"):
    """
    Run downloadCutout.py on the generated list.

    Tries:
      1) uv run downloadCutout.py --list {list_path}
      2) python3 downloadCutout.py --list {list_path}
      3) {sys.executable} downloadCutout.py --list {list_path}
    """
    print(f"[HSC] Launching downloader with list: {list_path}")

    # 1) uv run
    try:
        completed = subprocess.run(
            ["uv", "run", script_path, "--list", list_path],
            check=True
        )
        print("[HSC] Finished via: uv run")
        return
    except FileNotFoundError:
        print("[HSC] 'uv' not found; falling back to python3...")
    except subprocess.CalledProcessError as e:
        print(f"[HSC] 'uv run' failed (return code {e.returncode}); trying python3...")

    # 2) python3
    try:
        completed = subprocess.run(
            ["python3", script_path, "--list", list_path],
            check=True
        )
        print("[HSC] Finished via: python3")
        return
    except FileNotFoundError:
        print("[HSC] 'python3' not found; trying current interpreter...")
    except subprocess.CalledProcessError as e:
        print(f"[HSC] 'python3' failed (return code {e.returncode}); trying current interpreter...")

    # 3) current interpreter
    try:
        completed = subprocess.run(
            [sys.executable, script_path, "--list", list_path],
            check=True
        )
        print(f"[HSC] Finished via: {sys.executable}")
    except Exception as e:
        print(f"[HSC] All launch methods failed: {e}")
        raise


# -------------------------
# CLI
# -------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Build HSC list, optionally run downloadCutout.py, and optionally fetch radio cutouts.")
    p.add_argument("--survey", required=True, help="Survey name (e.g. 'vlass')")
    p.add_argument("--dropout", required=True, help="Dropout band (e.g. 'i')")
    p.add_argument("--root", default="/home/kongyw/hsc_wide_dropout_2510",
                   help="Root output directory (default: /home/kongyw/hsc_wide_dropout_2510)")
    p.add_argument("--catalog", help="Path to catalog CSV (default: {root}/{survey}/{survey}_{dropout}_dropout.csv)")
    p.add_argument("--filters", nargs="+", default=['g','r','i','z','y'], help="HSC filters")
    p.add_argument("--cutout-sizes", nargs="+", type=float, default=[2,4,8,12,30,90,150], help="HSC semi-sizes in arcsec")
    p.add_argument("--rerun", default="s23b_wide", help="HSC rerun (default s23b_wide)")

    p.add_argument("--radio-sizes", nargs="+", type=float, default=[1,3,5], help="Radio cutout sizes in arcmin (LoTSS/VLASS)")
    p.add_argument("--fetch-lotss", action="store_true", help="Fetch LoTSS cutouts")
    p.add_argument("--fetch-vlass", action="store_true", help="Fetch VLASS cutouts (requires astroquery)")

    p.add_argument("--download-script", default="downloadCutout.py",
                   help="Path to downloadCutout.py (default: downloadCutout.py)")
    p.add_argument("--run-download", action="store_true",
                   help="If set, run downloadCutout.py after writing the HSC list (default: do not run).")
    return p.parse_args()


def main():
    args = parse_args()

    survey = args.survey
    dropout = args.dropout
    root = os.path.abspath(args.root)
    os.makedirs(root, exist_ok=True)

    # default catalog UNDER root
    catalog_path = args.catalog or os.path.join(root, f"{survey}/{survey}_{dropout}_dropout.csv")
    if not os.path.isfile(catalog_path):
        print(f"[ERROR] Catalog not found: {catalog_path}")
        sys.exit(1)

    # 1) Build HSC list
    list_path = build_hsc_cutout_list(
        survey, dropout, root,
        catalog_path=catalog_path,
        cutout_sizes=tuple(args.cutout_sizes),
        filters=tuple(args.filters),
        rerun=args.rerun,
        output_txt=None,  # -> {root}/{survey}_{dropout}.txt
    )

    # 2) Optionally run downloadCutout.py
    if args.run_download:
        run_download_cutout(list_path, script_path=args.download_script)
    else:
        print("[HSC] Skipping downloadCutout.py run (use --run-download to enable).")

    # 3) Optionally fetch LoTSS (under root/{survey}/{survey}_{dropout})
    if args.fetch_lotss:
        lotss_out = os.path.join(root, f"{survey}/{survey}_{dropout}")
        os.makedirs(lotss_out, exist_ok=True)
        fetch_lotss_cutouts(
            survey, dropout,
            catalog_path=catalog_path,
            radio_sizes_arcmin=tuple(args.radio_sizes),
            output_dir=lotss_out,
        )

    # 4) Optionally fetch VLASS via CADC (under root/{survey}/{survey}_{dropout})
    if args.fetch_vlass:
        vlass_out = os.path.join(root, f"{survey}/{survey}_{dropout}")
        os.makedirs(vlass_out, exist_ok=True)
        fetch_vlass_cutouts_cadc(
            catalog_path=catalog_path,
            out_dir=vlass_out,
            radio_sizes_arcmin=tuple(args.radio_sizes),
            require_astroquery=True,
        )

    print("[DONE] Helper finished.")


if __name__ == "__main__":
    main()
