import os
import shutil
import subprocess
import urllib.request
import urllib.error
from pathlib import Path

from openmm.app import PDBFile, Modeller
from pdbfixer import PDBFixer


PROJECT_DIR = Path(__file__).resolve().parent

INPUT_FOLDER = PROJECT_DIR
FIXED_FOLDER = PROJECT_DIR
PDBQT_FOLDER = PROJECT_DIR

PREPARE_RECEPTOR = r"C:\Program Files (x86)\MGLTools-1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_receptor4.py"

MGL_PYTHON = r"C:\Program Files (x86)\MGLTools-1.5.7\python.exe"

PH = 7.4
ADD_MISSING_RESIDUES = False


def get_user_input():
    print()
    print("=" * 60)
    print("Protein Preparation")
    print("=" * 60)

    download = input("Download protein from database? (Y/N): ").strip().upper()
    if download == "Y":
        download_from_database = True
        while True:
            pdb_id = input("Enter PDB ID: ").strip().upper()
            if len(pdb_id) == 4 and pdb_id.isalnum():
                break
            print("Please enter a valid 4-character PDB ID.")

        redo_choice = input("Check PDB-REDO first? (Y/N): ").strip().upper()
        use_pdb_redo = redo_choice == "Y"

    else:
        download_from_database = False
        pdb_id = None
        use_pdb_redo = False

    chain_text = input("Enter chain IDs separated by commas (leave blank for all chains): ").strip()

    if chain_text == "":
        selected_chains = None

    else:
        selected_chains = [
            chain.strip().upper()
            for chain in chain_text.split(",")
            if chain.strip()
        ]

    return (download_from_database, pdb_id, use_pdb_redo, selected_chains)

def download_file(url, output_file):
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(request, timeout=60) as response:
            data = response.read()

        if not data:
            return False
        output_file.write_bytes(data)
        return True

    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError):
        if output_file.exists():
            output_file.unlink()
        return False

def validate_pdb_file(pdb_file):
    if not pdb_file.exists():
        return False

    if pdb_file.stat().st_size < 100:
        return False

    try:
        with open(pdb_file, "r", errors="ignore") as file:
            for line in file:
                if line.startswith(("ATOM", "HETATM", "HEADER", "MODEL")):
                    return True

    except OSError:
        return False

    return False

def download_from_pdb_redo(pdb_id):
    pdb_id = pdb_id.lower()
    output_file = PROJECT_DIR / f"{pdb_id.upper()}_PDBREDO.pdb"
    url = (
        f"https://pdb-redo.eu/db/"
        f"{pdb_id}/{pdb_id}_final.pdb"
    )
    print()
    print("Checking PDB-REDO...")
    print(f"URL: {url}")

    success = download_file(url, output_file)

    if not success:
        print("PDB-REDO PDB file is not available.")
        return None

    if not validate_pdb_file(output_file):
        print("Downloaded PDB-REDO file is not a valid PDB file.")

        if output_file.exists():
            output_file.unlink()
        return None
    print("PDB-REDO structure found.")
    print(f"Downloaded: {output_file.name}")
    return output_file

def download_from_rcsb(pdb_id):
    pdb_id = pdb_id.upper()
    output_file = PROJECT_DIR / f"{pdb_id}.pdb"
    if output_file.exists():
        print(f"{output_file} already exists.")
        return output_file
    url = (
        f"https://files.rcsb.org/download/"
        f"{pdb_id}.pdb"
    )
    print()
    print(f"Downloading {pdb_id} from RCSB PDB...")
    print(f"URL: {url}")
    success = download_file(url, output_file)

    if not success:
        raise RuntimeError(f"Could not download {pdb_id} from RCSB PDB.")

    if not validate_pdb_file(output_file):
        if output_file.exists():
            output_file.unlink()

        raise RuntimeError(f"Downloaded file for {pdb_id} is not a valid PDB file.")

    print("RCSB download complete.")
    return output_file

def obtain_protein(pdb_id, use_pdb_redo):
    if use_pdb_redo:
        redo_file = download_from_pdb_redo(pdb_id)
        if redo_file is not None:
            return redo_file
        print()
        print("PDB-REDO structure unavailable.Falling back to RCSB PDB.")
    return download_from_rcsb(pdb_id)

def keep_selected_chains(input_pdb, selected_chains):
    if selected_chains is None:
        print("Using all chains.")
        return input_pdb
    print()
    print("Selecting chains: " + ", ".join(selected_chains))
    input_pdb = Path(input_pdb)
    fixer = PDBFixer(filename=str(input_pdb))
    chains_to_delete = []
    available_chains = []

    for chain in fixer.topology.chains():
        available_chains.append(chain.id)
        if chain.id not in selected_chains:
            chains_to_delete.append(chain)
    print("Available chains: " + ", ".join(available_chains))

    missing_requested_chains = [
        chain
        for chain in selected_chains
        if chain not in available_chains
    ]

    if missing_requested_chains:
        raise RuntimeError("Requested chain(s) not found: "+ ", ".join(missing_requested_chains))
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.delete(chains_to_delete)
    output_file = input_pdb.with_name(input_pdb.stem + "_selected.pdb")
    with open(output_file, "w") as file:
        PDBFile.writeFile(modeller.topology, modeller.positions, file, keepIds=True)
    print(f"Selected-chain structure created: {output_file.name}")
    return output_file

def remove_water_from_pdb(pdb_file):
    pdb_file = Path(pdb_file)

    water_names = {"HOH", "WAT", "H2O"}
    cleaned_lines = []
    with open(pdb_file, "r", errors="ignore") as file:
        for line in file:
            if line.startswith(("ATOM", "HETATM")):
                residue = line[17:20].strip()
                if residue in water_names:
                    continue
            cleaned_lines.append(line)

    with open(pdb_file, "w") as file:
        file.writelines(cleaned_lines)

def repair_protein(input_file):
    input_file = Path(input_file)
    print()
    print(f"Repairing {input_file.name}")
    fixer = PDBFixer(filename=str(input_file))
    fixer.findMissingResidues()
    if not ADD_MISSING_RESIDUES:
        fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(PH)
    output_file = FIXED_FOLDER / (input_file.stem + "_fixed.pdb")

    with open(output_file, "w") as file:
        PDBFile.writeFile(fixer.topology, fixer.positions, file, keepIds=True)
    remove_water_from_pdb(output_file)
    print(f"Fixed protein created: {output_file.name}")
    print("Water molecules removed.")
    return output_file

def prepare_receptor(fixed_pdb):
    fixed_pdb = Path(fixed_pdb)
    receptor_name = fixed_pdb.stem
    output_pdbqt = PDBQT_FOLDER / (receptor_name + ".pdbqt")
    command = [MGL_PYTHON, PREPARE_RECEPTOR, "-r", str(fixed_pdb), "-o", str(output_pdbqt), "-A", "hydrogens", "-U", "nphs_lps_waters"]

    try:
        print()
        print(f"Preparing receptor: {receptor_name}")
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        if result.stdout:
            print(result.stdout)
        if output_pdbqt.exists():
            print(f"PDBQT created: {output_pdbqt.name}")
            return output_pdbqt

        print("PDBQT file was not created.")
        return None

    except FileNotFoundError as error:
        print()
        print("MGLTools could not be found.")
        print(error)
        return None

    except subprocess.CalledProcessError as error:
        print()
        print("Failed to prepare receptor.")
        if error.stdout:
            print(error.stdout)
        if error.stderr:
            print(error.stderr)
        return None


def process_protein(protein_file):
    try:
        fixed_pdb = repair_protein(protein_file)
        pdbqt = prepare_receptor(fixed_pdb)
        if pdbqt is None:
            return False
        return True

    except Exception as error:
        print()
        print(f"Error processing: {protein_file}")
        print(error)
        return False

def get_protein_files():
    proteins = []
    for file in INPUT_FOLDER.iterdir():
        if not file.is_file():
            continue
        if file.suffix.lower() != ".pdb":
            continue
        name = file.stem.lower()
        if "_selected" in name:
            continue
        if "_fixed" in name:
            continue
        if "_pdbredo" in name:
            continue
        proteins.append(file)
    return sorted(proteins)

def process_all_proteins():
    proteins = get_protein_files()
    if len(proteins) == 0:
        print()
        print("No original PDB files found.")
        return
    successful = 0
    failed = 0
    print()
    print(f"Found {len(proteins)} protein(s).")
    for protein in proteins:
        print()
        print(f"Processing {protein.name}")
        success = process_protein(protein)
        if success:
            successful += 1
        else:
            failed += 1

    print()
    print("Protein preparation finished.")
    print(f"Successful : {successful}")
    print(f"Failed     : {failed}")

def print_summary():
    print()
    print("=" * 60)
    print("Protein Preparation Completed")
    print("=" * 60)
    print(f"Project directory : {PROJECT_DIR}")
    print(f"Fixed proteins    : {FIXED_FOLDER}")
    print(f"Prepared PDBQT    : {PDBQT_FOLDER}")
    print("=" * 60)

def main():
    (download_from_database, pdb_id, use_pdb_redo, selected_chains) = get_user_input()
    if download_from_database:
        protein = obtain_protein(pdb_id, use_pdb_redo)
        if selected_chains is not None:
            protein = keep_selected_chains(protein, selected_chains)
        process_protein(protein)

    else:
        process_all_proteins()
    print_summary()

if __name__ == "__main__":
    main()