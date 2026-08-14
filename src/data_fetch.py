"""
Data fetching utilities for GLUT pathogenicity analysis.
"""

from pathlib import Path
import gzip
import requests
import tqdm


def get_alphamissense(output_path: str | Path) -> None:
    """Download and extract AlphaMissense proteome-wide predictions from Zenodo."""
    output_path = Path(output_path)
    archive = output_path.parent / "archive.tsv.gz"
    if output_path.exists():
        print(f"AlphaMissense data already exists at: {output_path}")
        return

    url = "https://zenodo.org/records/10813168/files/AlphaMissense_aa_substitutions.tsv.gz?download=1"
    print(f"Downloading AlphaMissense data from {url} ...")
    response = requests.get(url, stream=True)
    response.raise_for_status()
    total_size = int(response.headers.get("content-length", 0))

    tqdm_params = {
        "desc": "Downloading AlphaMissense",
        "total": total_size,
        "miniters": 1,
        "unit": "B",
        "unit_scale": True,
        "unit_divisor": 1024,
    }
    with tqdm.tqdm(**tqdm_params) as pb:
        with open(archive, "wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                handle.write(chunk)
                pb.update(len(chunk))

    print("Download complete, uncompressing...")
    with gzip.open(archive, "rb") as f_in, open(output_path, "wb") as f_out:
        f_out.writelines(f_in)
    archive.unlink(missing_ok=True)
    print(f"AlphaMissense data saved to: {output_path}")


def get_fasta(uniprot_id: str, output_path: str | Path) -> None:
    """Download protein FASTA sequence from EBI Proteins API."""
    output_path = Path(output_path)
    if output_path.exists():
        return
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}"
    response = requests.get(url, headers={"accept": "text/x-fasta"})
    response.raise_for_status()
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(response.text)


def get_alphafold_pdb(uniprot_id: str, output_path: str | Path) -> None:
    """Download AlphaFold2 predicted structure for a UniProt ID."""
    output_path = Path(output_path)
    if output_path.exists():
        return
    url = f"https://alphafold.com/api/prediction/{uniprot_id}"
    records = requests.get(url, headers={"accept": "application/json"}).json()
    pdb_url = records[0]["pdbUrl"]
    pdb_resp = requests.get(pdb_url)
    pdb_resp.raise_for_status()
    with open(output_path, "wb") as f:
        f.write(pdb_resp.content)


def get_protein_length(uniprot_id: str) -> int:
    """Retrieve canonical protein sequence length from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/search?query=accession:{uniprot_id}"
    records = requests.get(url, headers={"accept": "application/json"}).json()
    return records["results"][0]["sequence"]["length"]
