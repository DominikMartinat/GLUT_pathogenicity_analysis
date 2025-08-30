from pathlib import Path
import subprocess


def pymissense(am_file, pdb_file, uniprot_id, output_dir):
    am_file = Path(am_file)
    pdb_file = Path(pdb_file)
    output_dir = Path(output_dir)
    output = output_dir / f'{uniprot_id}-edit.pdb'

    #skip if prediction already exists
    if output.exists():
        print(f"AlphaMissense prediction for {uniprot_id} already exists")
        return
    
    command = f"pymissense --tsv {am_file} --pdbpath {pdb_file} {uniprot_id} {output_dir}"
    subprocess.check_output(command, shell=True, text=True)
    output_dir.joinpath(uniprot_id+'.pdf').unlink(missing_ok=True)

def deepTMHMM(application, fasta_file, uniprot_id, out_dir):
    out_dir = Path(out_dir)
    output = out_dir / f'{uniprot_id}.3line'

    #skip if prediction already exists
    if output.exists():
        print(f"DeepTMHMM prediction for {uniprot_id} already exists")
        return
    job = application.cli(args=f'--fasta {fasta_file}')
    job.save_files(output_dir=out_dir, path_filter='*.3line')
    out_dir.joinpath('predicted_topologies.3line').replace(output)

