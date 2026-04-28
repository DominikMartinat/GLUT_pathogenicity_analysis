from pathlib import Path
import shutil
import zipfile
import pandas as pd

def mole_zip2json(file_path, output):
    temp_dir = output.parent / 'temp'
    Path.mkdir(temp_dir)
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)
    json_file = temp_dir / 'json/data.json'
    json_file.replace(output)
    shutil.rmtree(temp_dir)

def unite_csvs(filelist, output, id_col = 1, colnames = None, return_df=False):
    united_df = pd.read_csv(filelist[0])
    cols_added = 0
    for f in filelist[1:]:
        pathogenicities = []
        df = pd.read_csv(f)
        for index, row in united_df.iterrows():
            identifier = row.iloc[id_col-1]
            pat = df[df.iloc[:,id_col-1]==identifier]['average_pathogenicity']
            pathogenicities.append(pat.iloc[0])
        cols_added += 1
        united_df.loc[:, str(cols_added)] = pathogenicities
    if colnames:
        united_df.columns = colnames
    united_df.to_csv(output, index=False)
    if return_df:
        return united_df
