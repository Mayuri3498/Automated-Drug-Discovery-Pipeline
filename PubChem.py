import re
import time
import csv
import urllib.request
import urllib.error
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy


PROJECT_DIR = Path(__file__).resolve().parent
SUMMARY_FILE = PROJECT_DIR / "summary.txt"
INPUT_FOLDER = PROJECT_DIR / "Molecules"
PDB_FOLDER = PROJECT_DIR / "Ligands_PDB"
PDBQT_FOLDER = PROJECT_DIR / "Ligands_PDBQT"

DOWNLOAD_SUCCESS_CSV = PROJECT_DIR / "successful_downloads.csv"
DOWNLOAD_FAILED_CSV = PROJECT_DIR / "failed_downloads.csv"
PREPARATION_SUCCESS_CSV = PROJECT_DIR / "successful_preparations.csv"
PREPARATION_FAILED_CSV = PROJECT_DIR / "failed_preparations.csv"

GENERATE_3D = True
OPTIMIZE = True
FORCE_FIELD = "MMFF"

INPUT_FOLDER.mkdir(exist_ok=True)
PDB_FOLDER.mkdir(exist_ok=True)
PDBQT_FOLDER.mkdir(exist_ok=True)

def extract_cids():
    with open(SUMMARY_FILE, "r", encoding="utf-8", errors="ignore") as file:
        text = file.read()
    cids = re.findall(r"CID:\s*(\d+)", text, re.IGNORECASE)
    return sorted(set(cids), key=int)

def download_sdf(cid):
    output_file = INPUT_FOLDER / f"{cid}.sdf"
    if output_file.exists() and output_file.stat().st_size > 0:
        print(f"Already exists: {output_file.name}")
        return True
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/"
        f"rest/pug/compound/cid/{cid}/SDF"
        "?record_type=2d"
    )
    print(f"Downloading CID: {cid}")

    try:
        request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(request, timeout=60) as response:
            data = response.read()
        if not data:
            raise RuntimeError("Empty response received.")
        if b"$$$$" not in data:
            raise RuntimeError("Downloaded file is not a valid SDF.")
        with open(output_file, "wb") as file:
            file.write(data)
        print(f"Downloaded successfully: {cid}.sdf")
        return True

    except urllib.error.HTTPError as error:
        print(f"Failed CID {cid}: HTTP {error.code} - {error.reason}")
        return False

    except urllib.error.URLError as error:
        print(f"Failed CID {cid}: Connection error - {error.reason}")
        return False

    except Exception as error:
        print(f"Failed CID {cid}: {error}")
        return False
    finally:
        time.sleep(0.25)

def write_download_csv(filename, cids, status):
    with open(filename, "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["Molecule Number", "CID", "Filename", "Status"])

        for cid in cids:
            writer.writerow([cid, cid,f"{cid}.sdf", status])
    print(f"Created: {filename}")

def download_all():
    if not SUMMARY_FILE.exists():
        print(f"Summary file not found: {SUMMARY_FILE}")
        return [], []
    cids = extract_cids()
    print(f"\nCIDs found: {len(cids)}")

    if not cids:
        print("No CIDs found.")
        return [], []
    successful_cids = []
    failed_cids = []
    print("\nStarting downloads...\n")

    for number, cid in enumerate(cids, start=1):
        print(f"[{number}/{len(cids)}]")
        success = download_sdf(cid)
        if success:
            successful_cids.append(cid)
        else:
            failed_cids.append(cid)

    write_download_csv(DOWNLOAD_SUCCESS_CSV, successful_cids, "Successful")
    write_download_csv(DOWNLOAD_FAILED_CSV, failed_cids, "Failed")
    print("\nDownload Summary")
    print(f"Total CIDs : {len(cids)}")
    print(f"Successful : {len(successful_cids)}")
    print(f"Failed     : {len(failed_cids)}")
    return successful_cids, failed_cids

def get_sdf_files():
    return sorted(INPUT_FOLDER.glob("*.sdf"), key=lambda x: int(x.stem) if x.stem.isdigit() else x.stem)

def read_sdf(sdf_file):
    supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
    mol = supplier[0]
    if mol is None:
        raise RuntimeError(f"Cannot read {sdf_file}")
    return mol

def add_hydrogens(mol):
    return Chem.AddHs(mol, addCoords=True)

def generate_3d(mol):
    if not GENERATE_3D:
        return mol
    result = AllChem.EmbedMolecule(mol, randomSeed=42)

    if result != 0:
        raise RuntimeError("RDKit could not generate a 3D conformer.")
    return mol

def optimize_geometry(mol):
    if not OPTIMIZE:
        return mol

    if FORCE_FIELD.upper() == "UFF":
        result = AllChem.UFFOptimizeMolecule(mol)
    else:
        result = AllChem.MMFFOptimizeMolecule(mol)

    if result == -1:
        raise RuntimeError(f"{FORCE_FIELD} geometry optimization failed.")
    return mol

def save_pdb(mol, filename):
    Chem.MolToPDBFile(mol, str(filename))

def prepare_ligand(mol):
    preparator = MoleculePreparation()
    setup_list = preparator.prepare(mol)
    if len(setup_list) == 0:
        raise RuntimeError("Meeko could not prepare the ligand.")
    return setup_list[0]

def save_pdbqt(setup, output_file):
    result = PDBQTWriterLegacy.write_string(setup)
    if isinstance(result, tuple):
        pdbqt_string = result[0]
    else:
        pdbqt_string = result

    with open(output_file, "w", encoding="utf-8") as file:
        file.write(pdbqt_string)

def process_ligand(sdf_file):
    molecule_number = sdf_file.stem
    try:
        print(f"\nPreparing molecule: {molecule_number}")
        mol = read_sdf(sdf_file)
        mol = add_hydrogens(mol)
        mol = generate_3d(mol)
        mol = optimize_geometry(mol)
        pdb_file = PDB_FOLDER / f"{molecule_number}.pdb"
        save_pdb(mol, pdb_file)
        setup = prepare_ligand(mol)
        pdbqt_file = PDBQT_FOLDER / f"{molecule_number}.pdbqt"

        save_pdbqt(setup, pdbqt_file)
        print(f"Successfully prepared molecule: {molecule_number}")
        return True

    except Exception as error:
        print(f"\nFailed preparation: {molecule_number}")
        print(error)
        return False

def write_preparation_csv(filename, molecule_numbers, status):
    with open(filename, "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["Molecule Number", "SDF Filename", "PDB Filename", "PDBQT Filename", "Status"])
        for molecule_number in molecule_numbers:
            writer.writerow([molecule_number, f"{molecule_number}.sdf", f"{molecule_number}.pdb", f"{molecule_number}.pdbqt", status])
    print(f"Created: {filename}")

def process_all_ligands():
    ligands = get_sdf_files()
    if len(ligands) == 0:
        print("No SDF files found.")
        return [], []
    successful = []
    failed = []

    print("\n")
    print("=" * 60)
    print("LIGAND PREPARATION STARTED")
    print("=" * 60)

    for number, ligand in enumerate(ligands, start=1):
        print(f"\n[{number}/{len(ligands)}]")
        success = process_ligand(ligand)
        molecule_number = ligand.stem
        if success:
            successful.append(molecule_number)
        else:
            failed.append(molecule_number)

    write_preparation_csv(PREPARATION_SUCCESS_CSV, successful, "Successful")
    write_preparation_csv(PREPARATION_FAILED_CSV, failed, "Failed")
    print("\n")
    print("=" * 60)
    print("PREPARATION SUMMARY")
    print("=" * 60)
    print(f"Total molecules : {len(ligands)}")
    print(f"Successful      : {len(successful)}")
    print(f"Failed          : {len(failed)}")
    print(f"\nPDB files       : {PDB_FOLDER}")
    print(f"PDBQT files     : {PDBQT_FOLDER}")
    return successful, failed

def main():
    print("=" * 60)
    print("PUBCHEM DOWNLOAD AND LIGAND PREPARATION")
    print("=" * 60)
    successful_downloads, failed_downloads = download_all()
    print("\n")
    process_all_ligands()

    print("\n")
    print("=" * 60)
    print("WORKFLOW COMPLETED")
    print("=" * 60)
    print("\nGenerated files:")

    print(f"1. Successful downloads : {DOWNLOAD_SUCCESS_CSV.name}")
    print(f"2. Failed downloads     : {DOWNLOAD_FAILED_CSV.name}")
    print(f"3. Successful preparation : {PREPARATION_SUCCESS_CSV.name}")
    print(f"4. Failed preparation     : {PREPARATION_FAILED_CSV.name}")

if __name__ == "__main__":
    main()