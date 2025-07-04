# hsc_s23bwide_dropout

## Explanation of the Files in the Google Drive

The dropout candidates are stored in `.csv` files named using the format:  
`{survey_name}_{dropout}_drop_candidates.csv`

- **`{survey_name}`**: `'first'`, `'lotss'`, `'vlass'`  
- **`{dropout}`**: `'g'`, `'r'`, `'i'`, `'z'`

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
   - Sizes: 4″×4″, 8″×8″, and 12″×12″  
   - Centered on the **optical counterpart**, with the **radio source** labeled  
   - All filters are scaled consistently with optimized contrast for inspecting bright pixels  
   - The 8″ and 12″ cutouts also show the cross-match region  
   - Crossmatch radii:
     - HSC-FIRST: 1.0″  
     - HSC-VLASS: 1.5″  
     - HSC-LoTSS: 2.0″  

2. **LoTSS Radio Cutouts**  
   - Size: 1 arcmin  
   - Centered on the **radio source** for all HSC-LoTSS dropout candidates  
   - Annotated with:
     - Radio center  
     - Matched optical counterpart  
     - Crossmatch region  
