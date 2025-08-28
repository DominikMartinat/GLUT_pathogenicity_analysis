import pandas as pd

def residues(pymissnense_pdb, out_file):
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