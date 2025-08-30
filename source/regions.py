import pandas as pd
import json

def all_residues(pdb_file,out_file):
    residue_labels = []
    residue_names = []
    for line in open(pdb_file):
            l = line.split()
            id = l[0]
            if id == 'ATOM':
                atom_type = l[2]
                if atom_type == 'CA':
                    #append AA residue position
                    residue_labels.append(l[5])
                    #append AA residue name
                    residue_names.append(l[3])
    pd.DataFrame(data={'residue_label':residue_labels,'residue_name':residue_names}).to_csv(out_file, index = False)

#parse deeptmhmm output file to get intracellular, extracellular and membrane residues
def membrane_residues(deeptmhmm_file, out_dir, identifier, ext_cellular_p='O', int_cellular_p='I', membrane_p='M'):
    intracellular = {'residue_label':[],'residue_name':[]}
    extracellular = {'residue_label':[],'residue_name':[]}
    membrane = {'residue_label':[],'residue_name':[]}
    with open(deeptmhmm_file, "r") as f:
        lines = f.read().splitlines()
    sequence = lines[1].strip()
    annotation = lines[2].strip()
    for i in range(len(sequence)):
        region = annotation[i]
        if region == 'I':
            intracellular['residue_label'].append(i+1)
            intracellular['residue_name'].append(sequence[i])
        elif region == 'O':
            extracellular['residue_label'].append(i+1)
            extracellular['residue_name'].append(sequence[i])
        elif region == 'M':
            membrane['residue_label'].append(i+1)
            membrane['residue_name'].append(sequence[i])
    intracellular_file = out_dir / f'{int_cellular_p}_{identifier}.csv'
    extracellular_file = out_dir / f'{ext_cellular_p}_{identifier}.csv'
    membrane_file = out_dir / f'{membrane_p}_{identifier}.csv'

    pd.DataFrame(data=intracellular).to_csv(intracellular_file,index=False)
    pd.DataFrame(data=extracellular).to_csv(extracellular_file,index=False)
    pd.DataFrame(data=membrane).to_csv(membrane_file,index=False)

def binding_pockets(prankweb_csv, out_file):
    df = pd.read_csv(prankweb_csv)
    df = df[df[' pocket']==1][[' residue_label', ' residue_name']]
    df.columns = ['residue_label', 'residue_name']
    df.sort_values(by = 'residue_label').to_csv(out_file, index=False)

def lining_residues(moleonline_json, out_file):
    with open(moleonline_json, "r", encoding="utf-8") as f:
        data = json.load(f) 
    lining_res_dict = dict()
    lining_res = data['Channels']['Paths'][0]['Layers']['ResidueFlow']
    for lr in lining_res:
        lr_l = lr.split()
        lining_res_dict[lr_l[1]] = lr_l[0]
    lining_res_df = pd.DataFrame(data = {'residue_label':lining_res_dict.keys(),'residue_name':lining_res_dict.values()})
    lining_res_df = lining_res_df.astype({'residue_label':int})
    lining_res_df.sort_values(by='residue_label').to_csv(out_file, index=False)

def nonbinding_pocket_lining_residues(bp_residues_csv, l_residues_csv,out_file):
    lr_df = pd.read_csv(l_residues_csv)
    bp_df = pd.read_csv(bp_residues_csv)

    residue_labels = []
    residue_names = []
    for index, row in lr_df.iterrows():
        if row['residue_label'] not in bp_df['residue_label'].values:
            residue_labels.append(row['residue_label'])
            residue_names.append(row['residue_name'])
    pd.DataFrame(data = {'residue_label':residue_labels, 'residue_name':residue_names}).to_csv(out_file, index=False)