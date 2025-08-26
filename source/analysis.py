from pathlib import Path
import pandas as pd

RES_DIR = Path.cwd().joinpath('results')


def alphamissense_avg_pat(uniprot_id, verbose=False):
    pat_sum = 0
    count = 0
    file = RES_DIR.joinpath('pymissense').joinpath(uniprot_id+'-edit.pdb')
    for line in open(f'{file}'):
            l = line.split()
            id = l[0]
            if id == 'ATOM':
                atom_type = l[2]
                if atom_type == 'CA':
                    pat_sum += float(l[10])
                    count += 1

    avg = pat_sum/count
    if verbose:
        print(f'Average aa pathogenicity for {uniprot_id}: {round(avg,2)}')
    return avg

def deepTMHMM3line_to_csv(input_file, output_file):
    with open(input_file, "r") as f:
        lines = f.read().splitlines()

    sequence = lines[1].strip()
    annotations = lines[2].strip()

    # 2. Control of the protein lenght
    if len(sequence) != len(annotations):
        raise ValueError("Leght of sequence and the annotation are not the same!")

    # 3. Creating a list I (in the cell), M (membrane region), O (out from the cell)
    I_list, M_list, O_list = [], [], []

    for i, (aa, tag) in enumerate(zip(sequence, annotations), start=1):
        entry = {"position": i, "name": aa}
        if tag == 'I':
            I_list.append(entry)
        elif tag == 'M':
            M_list.append(entry)
        elif tag == 'O':
            O_list.append(entry)

    # 4. DataFrames creation and alignment by number of lines
    max_len = max(len(I_list), len(M_list), len(O_list))

    def pad_list(lst):
        return lst + [{"position": "", "name": ""}] * (max_len - len(lst))

    I_df = pd.DataFrame(pad_list(I_list))
    M_df = pd.DataFrame(pad_list(M_list))
    O_df = pd.DataFrame(pad_list(O_list))

    # 5. Renaming columns. Position - position of the amino acid in the sequence, Name - type of the amino acid in alphabetical abbreviations
    I_df.columns = ["I_position", "I_name"]
    M_df.columns = ["M_position", "M_name"]
    O_df.columns = ["O_position", "O_name"]

    # 6. Combining into one table
    final_df = pd.concat([O_df, M_df, I_df], axis=1)

    # 7. Save to CSV
    final_df.to_csv(output_file, index=False)

    print(f"The resulting file was saved as: {output_file}")

#assign pathogenicity to protein regions - I (in the cell), M (membrane region), O (out from the cell)
def region_pathogenicity(region_file, pdb_file):
    # Loading an Excel file with aa region predictions
    df = pd.read_csv(region_file)

    df.columns = [col.strip() for col in df.columns]

    # Converting positions to integers (if they are floats or text)
    for col in ['O_position', 'M_position', 'I_position']:
        df[col] = pd.to_numeric(df[col],downcast='integer')

    # Loading a PDB file counted by PyMissense
    with open(pdb_file, "r") as f:
        pdb_lines = f.readlines()

    # Parsing residues and pathogenicity from column 11
    residue_patho = {}
    for line in pdb_lines:
        l = line.split()
        if l[0] == "ATOM":
                resnum = int(l[5])
                patho = float(l[10])
                #take only first atom per aa
                if resnum not in residue_patho:
                    residue_patho[resnum] = patho

    # Function to obtain pathogenicity for a given position
    position_patho = lambda pos : residue_patho.get(int(pos)) if pd.notnull(pos) else None

    # Add pathogenicity to all three positions
    df['O_pathogenicity'] = df['O_position'].apply(position_patho)
    df['M_pathogenicity'] = df['M_position'].apply(position_patho)
    df['I_pathogenicity'] = df['I_position'].apply(position_patho)

    # Output
    return df
