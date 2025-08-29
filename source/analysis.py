import pandas as pd

def assign_pathogenicity(patho_csv, region_csv, out_file):
    patho_df = pd.read_csv(patho_csv)
    region_df = pd.read_csv(region_csv)
    pathogenicities = []

    for index, row in patho_df.iterrows():
        if row['residue_label'] in region_df['residue_label'].values:
            pathogenicities.append(row['pathogenicity'])
    region_df.insert(2,'pathogenicity',pathogenicities)
    region_df.to_csv(out_file,index=False)

def average_pathogenicity(patho_dir, region_prefix, up_id_list, out_file):
    avg_pathogenicities = []
    up_ids = []
    for up_id in up_id_list:
        patho_file = patho_dir / f'{region_prefix}_{up_id}.csv'
        df = pd.read_csv(patho_file)
        avg_pathogenicities.append(df['pathogenicity'].mean())
        up_ids.append(up_id)
    pd.DataFrame(data={'uniprot_id':up_ids,'average_pathogenicity':avg_pathogenicities}).to_csv(out_file,index=False)
        