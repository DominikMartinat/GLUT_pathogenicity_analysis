from pathlib import Path
import subprocess

RES_DIR = Path.cwd().joinpath('results')
DATA_DIR = Path.cwd().joinpath('data')

def pymissense(uniprot_id, delete_pdf=True):
    output_dir = RES_DIR.joinpath('pymissense')
    output = output_dir.joinpath(uniprot_id+'-edit.pdb')

    #skip if prediction already exists
    if output.exists():
        print(f"AlphaMissense prediction for {uniprot_id} already exists")
        return
    
    am_file = DATA_DIR.joinpath('AlphaMissense_aa_substitutions.tsv')
    pdb_file = DATA_DIR.joinpath('pdb').joinpath(uniprot_id+'.pdb')
    command = f"pymissense --tsv {am_file} --pdbpath {pdb_file} {uniprot_id} {output_dir}"
    subprocess.check_output(command, shell=True, text=True)
    output_dir.joinpath(uniprot_id+'.pdf').unlink(missing_ok=True)

def deepTMHMM(uniprot_id, application):
    out_dir = RES_DIR.joinpath('deeptmhmm')
    output = out_dir.joinpath(uniprot_id+'.3line')

    #skip if prediction already exists
    if output.exists():
        print(f"DeepTMHMM prediction for {uniprot_id} already exists")
        return

    fasta = DATA_DIR.joinpath('fasta').joinpath(uniprot_id+'.fasta')
    job = application.cli(args=f'--fasta {fasta}')
    job.save_files(output_dir=out_dir, path_filter='*.3line')
    out_dir.joinpath('predicted_topologies.3line').replace(output)

