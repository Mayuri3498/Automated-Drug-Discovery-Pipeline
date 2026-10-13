import subprocess
import csv
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent

OUTPUT_FOLDER = PROJECT_DIR / "Docking_Results"

DOCKED_CSV = PROJECT_DIR / "docked_molecules.csv"
UNDOCKED_CSV = PROJECT_DIR / "undocked_molecules.csv"

print("=" * 60)
print("        AutoDock Vina Batch Docking")
print("=" * 60)

protein_name = input("\nEnter receptor PDBQT file name\n").strip()
ligand_folder_name = input("\nEnter folder containing ligand PDBQT files\n").strip()

config_name = input("\nEnter Vina configuration file"
    "\n(Press Enter for dock.txt): ").strip()

if config_name == "":
    config_name = "dock.txt"

PROTEIN_FILE = PROJECT_DIR / protein_name
LIGAND_FOLDER = PROJECT_DIR / ligand_folder_name
CONFIG_FILE = PROJECT_DIR / config_name
OUTPUT_FOLDER.mkdir(exist_ok=True)

if not PROTEIN_FILE.exists():
    raise FileNotFoundError("\nProtein receptor not found:\n{PROTEIN_FILE}")

if not LIGAND_FOLDER.exists():
    raise FileNotFoundError(f"\nLigand folder not found:\n{LIGAND_FOLDER}")

if not CONFIG_FILE.exists():
    raise FileNotFoundError(f"\nConfiguration file not found:\n{CONFIG_FILE}")

ligand_files = sorted(LIGAND_FOLDER.glob("*.pdbqt"))

if len(ligand_files) == 0:
    raise RuntimeError(f"\n No ligand PDBQT files found in \n {LIGAND_FOLDER}")

print("\nInput Summary")
print("-" * 60)
print(f"Protein      : {PROTEIN_FILE.name}")
print(f"Ligand Folder: {LIGAND_FOLDER.name}")
print(f"Ligands      : {len(ligand_files)}")
print(f"Config File  : {CONFIG_FILE.name}")
print(f"Output Folder: {OUTPUT_FOLDER}")
print("-" * 60)

def is_already_docked(ligand_name):
    ligand_output_folder = (OUTPUT_FOLDER / ligand_name)
    output_pdbqt = (ligand_output_folder / "out.pdbqt")
    if not output_pdbqt.exists():
        return False

    if output_pdbqt.stat().st_size == 0:
        return False

    return True

def dock_ligand(ligand_file):
    ligand_name = ligand_file.stem
    ligand_output_folder = (OUTPUT_FOLDER / ligand_name)
    ligand_output_folder.mkdir(exist_ok=True)
    output_pdbqt = (ligand_output_folder / "out.pdbqt")
    output_log = (ligand_output_folder / "log.txt")
    if is_already_docked(ligand_name):
        print(f"\nSkipping already docked molecule: {ligand_name}")
        print(f"Existing output: {output_pdbqt}")
        return "already_docked"

    print(f"\nProcessing ligand: {ligand_name}")
    command = ["vina", "--receptor", str(PROTEIN_FILE), "--ligand", str(ligand_file), "--config", str(CONFIG_FILE), "--out", str(output_pdbqt), "--log", str(output_log), "--cpu", "4"]

    try:
        subprocess.run(command, capture_output=True, text=True, check=True)
        if not output_pdbqt.exists():
            print(f"Docking failed: {ligand_name}")
            print("out.pdbqt was not generated.")
            return "failed"

        if output_pdbqt.stat().st_size == 0:
            print(f"Docking failed: {ligand_name}")
            print("Generated out.pdbqt is empty.")
            return "failed"

        print(f"Docking completed: {ligand_name}")
        print(f"Output: {output_pdbqt}")
        return "newly_docked"

    except subprocess.CalledProcessError as error:
        print(f"Docking failed: {ligand_name}")
        if error.stderr:
            print(error.stderr)
        return "failed"

def write_csv(filename, molecules):
    with open(filename, "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["Molecule_Number", "Molecule_File", "Status", "Output_File"])
        for molecule in molecules:
            writer.writerow([molecule["number"], molecule["file"], molecule["status"], molecule["output"]])

def dock_all_ligands():
    newly_docked = []
    already_docked = []
    failed = []
    print("=" * 60)
    for number, ligand in enumerate(ligand_files, start=1):
        print(f"\n[{number}/{len(ligand_files)}]")
        result = dock_ligand(ligand)
        ligand_name = ligand.stem
        output_file = (
            OUTPUT_FOLDER /
            ligand_name /
            "out.pdbqt"
        )

        molecule_data = {"number": ligand_name, "file": ligand.name, "status": "", "output": ""}
        if result == "newly_docked":
            molecule_data["status"] = ("Newly Docked")
            molecule_data["output"] = (str(output_file))
            newly_docked.append(molecule_data)

        elif result == "already_docked":
            molecule_data["status"] = ("Already Docked")
            molecule_data["output"] = (str(output_file))
            already_docked.append(molecule_data)

        else:
            molecule_data["status"] = ("Failed / Undocked")
            molecule_data["output"] = ""
            failed.append(molecule_data)

    all_docked = (newly_docked + already_docked)
    write_csv(DOCKED_CSV, all_docked)
    write_csv(UNDOCKED_CSV, failed)
    print("\n" + "=" * 60)
    print("Batch docking completed.")
    print("=" * 60)
    print(f"Newly docked : {len(newly_docked)}")
    print(f"Already docked/skipped: {len(already_docked)}")
    print(f"Failed/undocked: {len(failed)}")
    print(f"Total ligands: {len(ligand_files)}")
    print("\nCSV files created:")
    print(f"Docked molecules  : {DOCKED_CSV}")
    print(f"Undocked molecules: {UNDOCKED_CSV}")
    print(f"\nResults saved in: {OUTPUT_FOLDER}")

def check_vina():
    try:
        subprocess.run(["vina", "--help"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return True

    except Exception:
        return False

def print_summary():
    print("\n")
    print("=" * 60)
    print("Docking Summary")
    print("=" * 60)
    print(f"Protein File  : {PROTEIN_FILE.name}")
    print(f"Ligand Folder : {LIGAND_FOLDER}")
    print(f"Config File   : {CONFIG_FILE.name}")
    print(f"Output Folder : {OUTPUT_FOLDER}")
    print("=" * 60)

def main():
    if not check_vina():
        print("\nERROR:")
        print("AutoDock Vina executable was not found.")
        print("Please install AutoDock Vina and ensure 'vina' is available in your system PATH.")
        return

    dock_all_ligands()
    print_summary()

if __name__ == "__main__":
    main()