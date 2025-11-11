#!/usr/bin/env python3
# make_imaging.py
#
# Build PNG imaging panels for dropout candidates.
# - Inputs: --survey, --dropout (required)
# - By default: produce HSC-only imaging panels (top row with g,r,i,z,y)
# - If --include-radio is set: add a radio row (FIRST/LoTSS/VLASS/TGSS) when files exist
#
# Expected files on disk:
#   Catalog CSV: {root}/{survey}/{survey}_{dropout}_dropout.csv   (contains: tract, ra, dec, RA, DEC, optional *_cmodel_mag)
#   HSC FITS:    {root}/{survey}/{survey}_{dropout}_hsc/{tract}_{FILTER}_{SIZE}.fits   <-- SIZE is an INTEGER
#   Radio FITS:  {root}/{survey}/{survey}_{dropout}/tract_{tract}_{radio}_{radio_size}.fits  <-- radio_size ∈ {1,3,5} arcmin
#
# Output:
#   PNGs at {root}/{survey}/{survey}_{dropout}_hsc_images/tract_{tract}_hsc_{size}_all.png
#
# Notes:
# - “Adaptive” display by default: global sigma-clipped stretch using all HSC filters per target.
# - Radio row is optional (off by default). Enable with --include-radio.
# - Radio panels are soft-convolved to the listed survey beam size for consistent display.

import os
import sys
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ImageNormalize, SqrtStretch
from astropy.convolution import convolve, Gaussian2DKernel
from matplotlib.patches import Circle
from matplotlib.patheffects import withStroke
from scipy.ndimage import gaussian_filter

from photutils.aperture import CircularAperture, aperture_photometry

# -------------------------
# Defaults / constants
# -------------------------

DEFAULT_FILTERS = ['g', 'r', 'i', 'z', 'y']
DEFAULT_HSC_SIZES = [2, 4, 8, 12, 30, 90, 150]  # INTEGER semi-sizes in arcsec
CROSSMATCH_RADIUS_ARCSEC = 1.5
RADIO_BEAMS = {  # FWHM in arcsec for visual smoothing
    'tgss': 25.0,
    'vlass': 2.5,
    'first': 5.0,
    'lotss': 6.0,
}
RADIO_ALLOWED_SIZES = (1, 3, 5)  # arcmin

# -------------------------
# Small helpers
# -------------------------

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)
    return path

def wcs_pixel_scale_arcsec(header) -> float:
    """Return an approximate pixel scale [arcsec/pix] from a FITS header."""
    if 'CD1_1' in header and isinstance(header['CD1_1'], (float, int)):
        return abs(header['CD1_1']) * 3600.0
    if 'CDELT1' in header and isinstance(header['CDELT1'], (float, int)):
        return abs(header['CDELT1']) * 3600.0
    return 0.168  # safe fallback

def load_fits(path: str):
    data, header = fits.getdata(path, header=True)
    if data.ndim > 2:
        data = np.squeeze(data)
    return data, header

def add_center_marks(ax, x_opt, y_opt, x_rad, y_rad, pixel_scale, cutout_size_arcsec, crossmatch_radius_arcsec):
    ax.plot(x_opt, y_opt, marker='x', color='blue', markersize=8, mew=1.5, label='HSC centre')
    ax.plot(x_rad, y_rad, marker='+', color='red', markersize=8, mew=1.5, label='Radio centre')
    if cutout_size_arcsec != 2:
        radius_pix = crossmatch_radius_arcsec / pixel_scale
        circle = Circle((x_rad, y_rad), radius_pix, edgecolor='red', facecolor='none',
                        linewidth=1.2, linestyle='-')
        ax.add_patch(circle)

def compute_global_norm(hsc_paths):
    pixels = []
    for p in hsc_paths:
        if not os.path.exists(p):
            continue
        try:
            img, _ = load_fits(p)
            finite = img[np.isfinite(img)]
            if finite.size:
                pixels.append(finite)
        except Exception:
            continue

    if not pixels:
        return ImageNormalize(vmin=-1, vmax=5, stretch=SqrtStretch())

    allpix = np.concatenate(pixels)
    _, bkg_median, bkg_std = sigma_clipped_stats(allpix, sigma=3.0)
    vmin = bkg_median - 1.0 * bkg_std
    vmax = bkg_median + 5.0 * bkg_std
    return ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())

def format_mag_label(row, filt, snr=None):
    parts = []
    if snr is not None:
        parts.append(f"S/N={snr:.2f}")
    mag_col = f"{filt}_cmodel_mag"
    if mag_col in row and not pd.isna(row[mag_col]):
        parts.append(f"cmodel_mag={row[mag_col]:.2f}")
    else:
        parts.append("cmodel_mag>limit")
    return "\n".join(parts)

def choose_radio_size_arcmin(hsc_size_arcsec: int) -> int:
    """
    Map HSC semi-size (arcsec) -> radio cutout size (arcmin) in {1,3,5}
    Rule:
      - ≤ 30″  -> 1′
      - = 90″  -> 3′
      - = 150″ -> 5′
    """
    if hsc_size_arcsec <= 30:
        return 1
    if hsc_size_arcsec == 90:
        return 3
    if hsc_size_arcsec == 150:
        return 5
    # Fallback (shouldn't happen with provided defaults)
    return 1

# -------------------------
# Core imaging function
# -------------------------

def make_imaging(
    survey_name: str,
    dropout: str,
    root: str,
    *,
    include_radio: bool = False,
    radio_types=('tgss','vlass','first','lotss'),
    filters=DEFAULT_FILTERS,
    cutout_sizes=DEFAULT_HSC_SIZES,   # must be integers
):
    """
    Create imaging panels for each object in the catalog.
    HSC-only by default; if include_radio=True, add a radio row beneath HSC.

    - HSC FITS filenames: {root}/{survey}/{survey}_{dropout}_hsc/{tract}_{filter}_{size}.fits  (size is INT)
    - Radio FITS filenames: {root}/{survey}/{survey}_{dropout}/tract_{tract}_{radio}_{radio_size}.fits
                             where radio_size ∈ {1,3,5} (arcmin) chosen per HSC size.
    - Output PNGs: {root}/{survey}/{survey}_{dropout}_hsc_images/tract_{tract}_hsc_{size}_all.png
    """
    catalog_path = os.path.join(root, f"{survey_name}/{survey_name}_{dropout}_dropout.csv")
    if not os.path.isfile(catalog_path):
        print(f"[ERROR] Catalog not found: {catalog_path}")
        sys.exit(1)

    df = pd.read_csv(catalog_path)

    if 'tract' not in df.columns:
        print("[ERROR] 'tract' column is required in the catalog.")
        sys.exit(1)

    # enforce integer HSC sizes strictly
    cutout_sizes = [int(s) for s in cutout_sizes]

    ra_opt_col  = 'ra'  if 'ra'  in df.columns else 'RA'
    dec_opt_col = 'dec' if 'dec' in df.columns else 'DEC'
    ra_rad_col  = 'RA'  if 'RA'  in df.columns else ra_opt_col
    dec_rad_col = 'DEC' if 'DEC' in df.columns else dec_opt_col

    out_dir = ensure_dir(os.path.join(root, f"{survey_name}/{survey_name}_{dropout}_hsc_images"))
    hsc_dir = os.path.join(root, f"{survey_name}/{survey_name}_{dropout}_hsc")
    radio_dir = os.path.join(root, f"{survey_name}/{survey_name}_{dropout}")

    tract_counter = defaultdict(int)

    for _, row in df.iterrows():
        tract_base = str(row['tract'])
        n = tract_counter[tract_base]
        tract_counter[tract_base] += 1
        tract = tract_base if n == 0 else f"{tract_base}_{n}"

        ra_optical, dec_optical = float(row[ra_opt_col]), float(row[dec_opt_col])
        ra_radio,   dec_radio   = float(row[ra_rad_col]), float(row[dec_rad_col])

        coord_optical = SkyCoord(ra=ra_optical*u.deg, dec=dec_optical*u.deg)
        coord_radio   = SkyCoord(ra=ra_radio  *u.deg, dec=dec_radio  *u.deg)

        offset = coord_radio.separation(coord_optical)
        print(f"[{tract}] Positional offset: {offset.arcsec:.3f} arcsec")

        for size in cutout_sizes:  # size is INT
            # HSC files for normalization
            hsc_paths = [os.path.join(hsc_dir, f"{tract}_{filt}_{size}.fits") for filt in filters]
            norm = compute_global_norm(hsc_paths)

            # layout
            num_hsc = len(filters)
            if include_radio:
                num_rows = 2
                num_cols = max(num_hsc, len(radio_types))
            else:
                num_rows = 1
                num_cols = num_hsc

            figsize = (4 * num_cols, 8 if include_radio else 4.5)
            fig = plt.figure(figsize=figsize)
            panel_hsc = 1
            panel_radio = num_cols + 1

            # ---- HSC row ----
            for filt in filters:
                image_path = os.path.join(hsc_dir, f"{tract}_{filt}_{size}.fits")  # strictly int size
                if not os.path.exists(image_path):
                    print(f("[HSC] Missing: {image_path}"))
                    ax = fig.add_subplot(num_rows, num_cols, panel_hsc)
                    ax.set_axis_off()
                    panel_hsc += 1
                    continue

                try:
                    image_data, header = load_fits(image_path)
                    wcs = WCS(header)

                    pixscale = wcs_pixel_scale_arcsec(header)
                    x_opt, y_opt = wcs.world_to_pixel(coord_optical)
                    x_rad, y_rad = wcs.world_to_pixel(coord_radio)

                    if size in (2, 4):
                        image_data = gaussian_filter(image_data, sigma=1.0)

                    ax = fig.add_subplot(num_rows, num_cols, panel_hsc, projection=wcs)
                    panel_hsc += 1

                    ax.set_title(f"{filt}", fontsize=14)
                    ax.set_aspect('auto')
                    ax.set_xlabel('RA')
                    ax.set_ylabel('Dec')
                    ax.imshow(image_data, origin='lower', cmap='binary', norm=norm)

                    snr = None
                    if size == 4:
                        r_pix = 2.0 / pixscale
                        ap = CircularAperture((x_opt, y_opt), r=r_pix)
                        phot = aperture_photometry(image_data, ap)
                        finite = image_data[np.isfinite(image_data)]
                        if finite.size:
                            _, _, bkg_std = sigma_clipped_stats(finite, sigma=3.0)
                            flux = phot['aperture_sum'][0]
                            snr = flux / (bkg_std * np.sqrt(np.pi * (r_pix**2)))
                        ap.plot(ax, color='cyan', lw=1.3)

                    add_center_marks(ax, x_opt, y_opt, x_rad, y_rad, pixscale, size, CROSSMATCH_RADIUS_ARCSEC)

                    label_text = format_mag_label(row, filt, snr=snr)
                    ax.text(0.05, 0.05, label_text.strip(), transform=ax.transAxes,
                            color='orange', fontsize=12,
                            path_effects=[withStroke(linewidth=1.5, foreground='black')])

                except Exception as e:
                    print(f"[HSC] Error loading {image_path}: {e}")
                    ax = fig.add_subplot(num_rows, num_cols, panel_hsc)
                    ax.set_axis_off()
                    panel_hsc += 1
                    continue

            # ---- Radio row (optional) ----
            if include_radio:
                radio_size_arcmin = choose_radio_size_arcmin(size)  # pick 1, 3, or 5 (arcmin)
                if radio_size_arcmin not in RADIO_ALLOWED_SIZES:
                    radio_size_arcmin = 1

                for rname in radio_types:
                    image_path = os.path.join(radio_dir, f"tract_{tract}_{rname}_{radio_size_arcmin}.fits")
                    if not os.path.exists(image_path):
                        print(f"[RADIO] Missing: {image_path}")
                        ax = fig.add_subplot(num_rows, num_cols, panel_radio)
                        ax.set_axis_off()
                        panel_radio += 1
                        continue

                    try:
                        image_data, header = load_fits(image_path)
                        wcs_r = WCS(header, naxis=2)

                        x_opt, y_opt = wcs_r.world_to_pixel(coord_optical)
                        x_rad, y_rad = wcs_r.world_to_pixel(coord_radio)

                        pixscale = wcs_pixel_scale_arcsec(header)

                        # Beam-matched smoothing (visual)
                        fwhm_arcsec = RADIO_BEAMS.get(rname, 6.0)
                        sigma_pix = max(fwhm_arcsec / (2.355 * pixscale), 0.0)
                        if sigma_pix > 0:
                            kernel = Gaussian2DKernel(x_stddev=sigma_pix)
                            image_data = convolve(image_data, kernel, preserve_nan=True)

                        finite = image_data[np.isfinite(image_data)]
                        if finite.size:
                            _, bkg_med, bkg_std = sigma_clipped_stats(finite, sigma=3.0)
                            image_data = image_data - bkg_med
                            vmin, vmax = (1.0 * bkg_std, 5.0 * bkg_std)
                            norm_r = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())
                        else:
                            norm_r = ImageNormalize(vmin=-1, vmax=5, stretch=SqrtStretch())

                        ax = fig.add_subplot(num_rows, num_cols, panel_radio, projection=wcs_r)
                        panel_radio += 1

                        ax.set_title(f"{rname.upper()} {radio_size_arcmin}'", fontsize=12)
                        ax.set_aspect('auto')
                        ax.set_xlabel('RA')
                        ax.set_ylabel('Dec')
                        ax.imshow(image_data, origin='lower', cmap='gray', norm=norm_r)

                        r_pix = (fwhm_arcsec / pixscale)
                        ap = CircularAperture((x_rad, y_rad), r=r_pix)
                        flux = aperture_photometry(image_data, ap)['aperture_sum'][0] if finite.size else np.nan
                        snr = (flux / (bkg_std * np.sqrt(np.pi * (r_pix**2)))) if finite.size else np.nan

                        ax.plot(x_opt, y_opt, marker='x', color='blue', markersize=6)
                        ax.plot(x_rad, y_rad, marker='+', color='red', markersize=6)
                        ap.plot(ax, color='cyan', lw=1.2)

                        label = f"S/N={snr:.1f}" if np.isfinite(snr) else "S/N=nan"
                        if rname == survey_name and 'Total_flux' in row and pd.notna(row['Total_flux']):
                            label += f"\n{row['Total_flux']:.2f} mJy"
                        label += f"\nOffset={offset.arcsec:.3f}\""

                        ax.text(0.02, 0.05, label, transform=ax.transAxes, color='orange',
                                fontsize=12, path_effects=[withStroke(linewidth=2, foreground='black')])

                    except Exception as e:
                        print(f"[RADIO] Error loading {image_path}: {e}")
                        ax = fig.add_subplot(num_rows, num_cols, panel_radio)
                        ax.set_axis_off()
                        panel_radio += 1
                        continue

            fig.suptitle(
                f"Tract {tract} HSC-{survey_name.upper()} {dropout}-dropout  "
                f" HSC{' + Radio' if include_radio else ''}  (HSC size={size*2}\")",
                fontsize=16
            )
            plt.tight_layout()
            plt.subplots_adjust(top=0.88)

            out_png = os.path.join(out_dir, f"tract_{tract}_hsc_{size}_all.png")  # size INT in name
            plt.savefig(out_png, dpi=150)
            plt.close(fig)
            print(f"[SAVE] {out_png}")

    print("[DONE] Imaging complete.")

# -------------------------
# CLI
# -------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Create HSC imaging panels (and optionally radio) for a survey dropout.")
    p.add_argument("--survey", dest="survey_name", required=True, help="Survey name (e.g. 'vlass')")
    p.add_argument("--dropout", required=True, help="Dropout band (e.g. 'i')")
    p.add_argument("--root", default="/home/kongyw/hsc_wide_dropout_2510", help="Root path containing survey folders (default: current dir)")

    # HSC rendering options — sizes are strictly INTEGER
    p.add_argument("--filters", nargs="+", default=DEFAULT_FILTERS, help="HSC filters to render")
    p.add_argument("--sizes", nargs="+", type=int, default=DEFAULT_HSC_SIZES, help="HSC semi-sizes (arcsec, INTEGER) used in filenames")

    # Radio (optional)
    p.add_argument("--include-radio", action="store_true", help="Include a radio row under HSC row")
    p.add_argument("--radio-types", nargs="+", default=['tgss','vlass','first','lotss'],
                   help="Radio surveys to show if files exist (radio size auto-chosen per HSC size)")
    return p.parse_args()

def main():
    args = parse_args()
    make_imaging(
        survey_name=args.survey_name,
        dropout=args.dropout,
        root=os.path.abspath(args.root),
        include_radio=args.include_radio,
        radio_types=tuple(args.radio_types),
        filters=tuple(args.filters),
        cutout_sizes=list(args.sizes),  # already ints
    )

if __name__ == "__main__":
    main()
