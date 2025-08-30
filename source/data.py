import requests
import os
import gzip
from pathlib import Path

DATA_DIR = Path.cwd().joinpath('data')

def get_alphamissense(output):
    #output_name = "AlphaMissense_aa_substitutions.tsv"
    #output = DATA_DIR.joinpath(output_name)
    output = Path(output)
    archive = output.parent / 'archive.gz'
    if output.exists():
        print("AlphaMissense data already exists")
        return
    request = requests.get("https://zenodo.org/records/10813168/files/AlphaMissense_aa_substitutions.tsv.gz?download=1", stream=True)
    with open(archive, "wb") as handle:
        for data in request.iter_content(chunk_size=10 * 1024):
            handle.write(data)
    unziped = gzip.open(archive, 'rb')
    with open(output, 'wb') as f_out:
        f_out.writelines(unziped)
    unziped.close()
    os.remove(archive)
    print ('AlphaMissense data downloaded')

def get_pdb(up_id, output):
    #pdb_dir = DATA_DIR.joinpath('pdb')
    #output = pdb_dir.joinpath(up_id+'.pdb')
    if os.path.exists(output):
        return
    records = requests.get(f"https://alphafold.com/api/prediction/{up_id}", headers={"accept": "application/json"}).json()
    pdbUrl = records[0]["pdbUrl"]
    with open(output, 'wb') as f:
        f.write(requests.get(pdbUrl).content)

def get_fasta(up_id, output):
    #fasta_dir = DATA_DIR.joinpath('fasta')
    #output = fasta_dir.joinpath(up_id+'.fasta')
    if os.path.exists(output):
        return
    records = requests.get(f"https://www.ebi.ac.uk/proteins/api/proteins/{up_id}", headers={"accept": "text/x-fasta"})
    with open(output, 'w') as f:
        f.write(records.text)