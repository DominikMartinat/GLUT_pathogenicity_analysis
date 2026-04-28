import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def assign_pathogenicity(patho_csv, region_csv, out_file):
    patho_df = pd.read_csv(patho_csv)
    region_df = pd.read_csv(region_csv)
    pathogenicities = []
    residue_labels = []
    residue_names = []

    for index, row in patho_df.iterrows():
        if row['residue_label'] in region_df['residue_label'].values:
            pathogenicities.append(row['pathogenicity'])
            residue_labels.append(row['residue_label'])
            residue_names.append(row['residue_name'])

    pd.DataFrame(data={'residue_label':residue_labels,'residue_name':residue_names, 'pathogenicity':pathogenicities}).to_csv(out_file,index=False)

def average_pathogenicity(patho_dir, region_prefix, up_id_list, out_file):
    avg_pathogenicities = []
    up_ids = []
    for up_id in up_id_list:
        patho_file = patho_dir / f'{region_prefix}_{up_id}.csv'
        df = pd.read_csv(patho_file)
        avg_pathogenicities.append(df['pathogenicity'].mean())
        up_ids.append(up_id)
    pd.DataFrame(data={'uniprot_id':up_ids,'average_pathogenicity':avg_pathogenicities}).to_csv(out_file,index=False)

def crate_heatmap(input_file, method_name, output_img, inverse_color=False):
    # Loading file
    df = pd.read_csv(input_file)

    # Setting 'filename' as index
    df.set_index("identifier", inplace=True)

    # Creating a heat map with values
    plt.figure(figsize=(10, len(df) * 0.4))

    if inverse_color:
        sns.heatmap(
            df,
            cmap="coolwarm_r",
            vmin=0,
            #vmax=1,
            linewidths=0.5,
            linecolor='gray',
            annot=True,
            fmt=".3f"
        )

    else:
        sns.heatmap(
        df,
        cmap="coolwarm",
        vmin=0,
        vmax=1,
        linewidths=0.5,
        linecolor='gray',
        annot=True,
        fmt=".3f"  
        )


    plt.title(f"GLUTs pathogenicity profile {method_name}", fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()


    plt.savefig(output_img, dpi=300, bbox_inches='tight')

    plt.show()