import os
import numpy as np
from astropy.table import Table

def chunk_array(arr, chunk_size):
    for i in range(0, len(arr), chunk_size):
        yield arr[i:i + chunk_size]

def main():
    # Get user input
    input_file = input("Enter the path to the input FITS file: ").strip()
    output_dir = input("Enter the path to the output directory: ").strip()

    # Validate input file
    if not os.path.isfile(input_file):
        print(f"❌ File not found: {input_file}")
        return

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Read input data
    df = Table.read(input_file).to_pandas()
    db_ids = df['db_id']

    # SQL templates
    sql_prefix = """SELECT
         main.object_id
        ,main.ra
        ,main.dec
        ,main.tract
        ,main.patch
        ,main.g_cmodel_mag
        ,main.g_cmodel_magerr
        ,main.g_cmodel_flux
        ,main.g_cmodel_fluxerr
        ,main.r_cmodel_mag
        ,main.r_cmodel_magerr
        ,main.r_cmodel_flux
        ,main.r_cmodel_fluxerr
        ,main.i_cmodel_mag
        ,main.i_cmodel_magerr
        ,main.i_cmodel_flux
        ,main.i_cmodel_fluxerr
        ,main.z_cmodel_mag
        ,main.z_cmodel_magerr
        ,main.z_cmodel_flux
        ,main.z_cmodel_fluxerr
        ,main.y_cmodel_mag
        ,main.y_cmodel_magerr
        ,main.y_cmodel_flux
        ,main.y_cmodel_fluxerr
        ,forced4.g_convolvedflux_0_20_mag
        ,forced4.g_convolvedflux_0_20_magerr
        ,forced4.r_convolvedflux_0_20_mag
        ,forced4.r_convolvedflux_0_20_magerr
        ,forced4.i_convolvedflux_0_20_mag
        ,forced4.i_convolvedflux_0_20_magerr
        ,forced4.z_convolvedflux_0_20_mag
        ,forced4.z_convolvedflux_0_20_magerr
        ,forced4.y_convolvedflux_0_20_mag
        ,forced4.y_convolvedflux_0_20_magerr
        ,forced4.g_convolvedflux_0_20_flag
        ,forced4.r_convolvedflux_0_20_flag
        ,forced4.i_convolvedflux_0_20_flag
        ,forced4.z_convolvedflux_0_20_flag
        ,forced4.y_convolvedflux_0_20_flag
        ,forced4.g_convolvedflux_0_20_flux
        ,forced4.g_convolvedflux_0_20_fluxerr
        ,forced4.r_convolvedflux_0_20_flux
        ,forced4.r_convolvedflux_0_20_fluxerr
        ,forced4.i_convolvedflux_0_20_flux
        ,forced4.i_convolvedflux_0_20_fluxerr
        ,forced4.z_convolvedflux_0_20_flux
        ,forced4.z_convolvedflux_0_20_fluxerr
        ,forced4.y_convolvedflux_0_20_flux
        ,forced4.y_convolvedflux_0_20_fluxerr
        ,main.a_g
        ,main.a_r
        ,main.a_i
        ,main.a_z
        ,main.a_y
        ,main.g_cmodel_flag, main.r_cmodel_flag, main.i_cmodel_flag, main.z_cmodel_flag, main.y_cmodel_flag
        ,main.merge_peak_g,main.merge_peak_r,main.merge_peak_i, main.merge_peak_z, main.merge_peak_y, 
        meas2.r_sdsscentroid_flag, meas2.i_sdsscentroid_flag, meas2.z_sdsscentroid_flag, meas2.y_sdsscentroid_flag,
        masks.g_mask_brightstar_any, masks.r_mask_brightstar_any, masks.i_mask_brightstar_any, masks.z_mask_brightstar_any, masks.y_mask_brightstar_any,
        masks.g_mask_brightstar_ghost, masks.r_mask_brightstar_ghost, masks.i_mask_brightstar_ghost, masks.z_mask_brightstar_ghost, masks.y_mask_brightstar_ghost,
        masks.g_mask_brightstar_ghost15, masks.r_mask_brightstar_ghost15, masks.i_mask_brightstar_ghost15, masks.z_mask_brightstar_ghost15, masks.y_mask_brightstar_ghost15,
        masks.g_mask_brightstar_halo, masks.r_mask_brightstar_halo, masks.i_mask_brightstar_halo, masks.z_mask_brightstar_halo, masks.y_mask_brightstar_halo,
        masks.g_mask_brightstar_blooming, masks.r_mask_brightstar_blooming, masks.i_mask_brightstar_blooming, masks.z_mask_brightstar_blooming, masks.y_mask_brightstar_blooming, masks.y_mask_brightstar_channel_stop,
        meas2.r_blendedness_abs, meas2.i_blendedness_abs, meas2.z_blendedness_abs, meas2.y_blendedness_abs,
        main.g_pixelflags_bright_objectcenter, main.r_pixelflags_bright_objectcenter, main.i_pixelflags_bright_objectcenter, main.z_pixelflags_bright_objectcenter, main.y_pixelflags_bright_objectcenter
FROM
        s23b_wide.forced as main
        INNER JOIN s23b_wide.meas2 AS meas2 ON (main.object_id = meas2.object_id)
        INNER JOIN s23b_wide.forced4 as forced4 ON (main.object_id = forced4.object_id)
        INNER JOIN s23b_wide.masks AS masks ON (main.object_id = masks.object_id)
WHERE
main.object_id in ("""

    sql_suffix = """)
    AND main.isprimary = 't'
    AND main.g_pixelflags_edge = 'f'
    AND main.r_pixelflags_edge = 'f'
    AND main.i_pixelflags_edge = 'f'
    AND main.z_pixelflags_edge = 'f'
    AND main.y_pixelflags_edge = 'f'
    AND main.g_pixelflags_saturatedcenter = 'f'
    AND main.r_pixelflags_saturatedcenter = 'f'
    AND main.i_pixelflags_saturatedcenter = 'f'
    AND main.z_pixelflags_saturatedcenter = 'f'
    AND main.y_pixelflags_saturatedcenter = 'f'
    AND main.g_pixelflags_bad = 'f'
    AND main.r_pixelflags_bad = 'f'
    AND main.i_pixelflags_bad = 'f'
    AND main.z_pixelflags_bad = 'f'
    AND main.y_pixelflags_bad = 'f'
    AND main.g_inputcount_value >= 3
    AND main.r_inputcount_value >= 3
    AND main.i_inputcount_value >= 5
    AND main.z_inputcount_value >= 5
    AND main.y_inputcount_value >= 5"""

    # Write each query
    for i, chunk in enumerate(chunk_array(db_ids, 100000), 1):
        chunk_formatted = ",\n".join(f"{id}" for id in chunk)
        full_query = f"-- Batch {i} {input_file} ---- \n" + sql_prefix + chunk_formatted + sql_suffix
        output_file = os.path.join(output_dir, f"{i}.sql")
        with open(output_file, "w") as f:
            f.write(full_query)

    print(f"✅ Generated {i} SQL files in {output_dir}/")

if __name__ == "__main__":
    main()
