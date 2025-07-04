# hsc_s23bwide_dropout
## Jupyter Notebook For Data Fetching
The dropout candidates are stored within the .csv files with the naming scheme `{survey_name}_{dropout}_drop_candidates.csv`.

The available inputs for `{survey_name}` are: 'first', 'lotss', and 'vlass'.

The available inputs for `{dropout}` are: 'g', 'r', 'i', 'z'.

The cutout files in the Google Drive are named by the 'tract' column in the data table, where 'tract' corresponds to its HSC position. However, it is quite common having multiple candidates within the same 'tract' as a single tract is ~1.5deg wide. For simplicity, I perform the following naming scheme just for simple distinguishment:

```python
tract_counter = defaultdict(int)

for index, row in data.iterrows():
    tract_base = str(row['tract'])
    count = tract_counter[tract_base]
    tract_counter[tract_base] += 1
    tract = tract_base if count == 0 else f"{tract_base}_{count}"

The jupyter notebook can be used blindly by just giving the survey name and the dropout input for fetching the available raw cutout fits files or processed images with the optical and radio centres labelled.
