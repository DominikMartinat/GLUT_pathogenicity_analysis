from pathlib import Path
import shutil
import zipfile

def mole_zip2json(file_path, output):
    temp_dir = output.parent / 'temp'
    Path.mkdir(temp_dir)
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)
    json_file = temp_dir / 'json/data.json'
    json_file.replace(output)
    shutil.rmtree(temp_dir)
