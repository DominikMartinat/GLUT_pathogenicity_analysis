import pandas as pd

def alphamissense_residues(pymissnense_pdb, out_file):
    residue_labels = []
    residue_names = []
    pathogenicities = []
    for line in open(pymissnense_pdb):
            l = line.split()
            id = l[0]
            if id == 'ATOM':
                atom_type = l[2]
                if atom_type == 'CA':
                    #append AA residue position
                    residue_labels.append(l[5])
                    #append AA residue name
                    residue_names.append(l[3])
                    #pathogenicity of AA
                    pathogenicities.append(l[10])
    pd.DataFrame(data={'residue_label':residue_labels,'residue_name':residue_names, 'pathogenicity':pathogenicities}).to_csv(out_file, index = False)

def polyphen2_residues(polyphen2_file, out_file):
    columns = [
        'ID', 'residue_label', 'residue_name', 'AA_alt', 'Training', 'Score', 'Prob', 'Class',
        'pathogenicity', 'PolyPhen_class', 'EVmutation_score', 'EVmutation_class'
    ]

    df = pd.read_csv(polyphen2_file, sep='\s+', names=columns, comment='#', na_values='nan')

    df = df[['residue_label', 'residue_name', 'pathogenicity']].copy()
    df = df.dropna(subset=['pathogenicity'])

    result = df.groupby(['residue_label', 'residue_name']).agg({'pathogenicity': 'mean'}).reset_index()
    result['pathogenicity'] = result['pathogenicity'].round(3)

    result.to_csv(out_file, index=False)

def sift_residues(sift_file, out_file):
    names = 'ACDEFGHIKLMNPQRSTVWY'
    active = False
    residue_labels = []
    residue_names = []
    pathogenicities = []
    
    with open(sift_file,'r') as f:
        for line in f:
            if active and line.startswith('pos'):
                continue
            elif active:
                l = line.split()
                position = l[0][:-1]
                name = l[0][-1]
                pathos = l[1:]
                c = 0
                for i in range(20):
                    if names[i] != name:
                        c += float(pathos[i])
                avg_pat = c/19
                residue_labels.append(position)
                residue_names.append(name)
                pathogenicities.append(avg_pat)
            elif line.startswith('pos'):
                active = True
    pd.DataFrame(data={'residue_label':residue_labels, 'residue_name':residue_names, 'pathogenicity':pathogenicities}).to_csv(out_file)