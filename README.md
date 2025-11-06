# hsc_s23bwide_dropout

## Explanation of the Helper scripts
The dropout candidates are stored in `.csv` files named using the format:  
`{survey_name}_{dropout}_drop_candidates.csv`

- **`{survey_name}`**: `'first'`, `'lotss'`, `'vlass'`  
- **`{dropout}`**: `'g'`, `'r'`, `'i'`, `'z'`

Using the fetch_all.py helper, one can obtain the raw '.fits' cutouts from HSC S23B Wide, LoTSS DR2, and the latest VLASS Epoch.
example input
run fetch_all.py --survey {survey_name} --dropout {dropout}
optional for getting the cutouts: --fetch-lotss --fetch-vlass --run-download

The cutout files are named based on the `tract` column in the data table, which corresponds to the HSC sky region. Since a single tract covers ~1.5 degrees, multiple candidates often fall within the same tract. To distinguish them, the following naming scheme is used:

```python
tract_counter = defaultdict(int)

for index, row in data.iterrows():
    tract_base = str(row['tract'])
    count = tract_counter[tract_base]
    tract_counter[tract_base] += 1
    tract = tract_base if count == 0 else f"{tract_base}_{count}"
```

### Available Cutouts

1. **HSC S23B Optical Cutouts**  
   - Sizes[arcsec]: [4, 8, 16, 24, 60, 180, 300]  
   - Centred on the **optical counterpart**, with the **radio source** labeled  
   - All filters are scaled consistently with optimized contrast for inspecting bright pixels  
   - The 8″ and 12″ cutouts also show the cross-match region  
   - Crossmatch radii:
     - HSC-FIRST: 1.0″  
     - HSC-VLASS: 1.5″  
     - HSC-LoTSS: 2.0″  

2. **VLASS & FIRST**
   - Sizes[arcmin]: [1, 3, 5]
   - Centred on the **radio source** for all dropout candidates
   - Annotated with:
     - Radio centre  
     - Matched optical counterpart  
     - Crossmatch region
3. **LoTSS DR2**
   - Matching cutout images for the HSC-VLASS dropout candidates

## SQL generating scripts
Applying the least selection flags to obtain the inital samples from the server-side crossmatch results between S23B wide and radio source catalogs.
