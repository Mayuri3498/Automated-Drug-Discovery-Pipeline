import os
import csv
import shutil

DOCKING_FOLDER = "Docking_Results"
OUTPUT_CSV = "binding_energies.csv"
CUTOFF_CSV = "cutoff_binding_energies.csv"
CUTOFF_FOLDER = "Cutoff_Results"

cutoff = float(input("Enter binding energy cutoff (kcal/mol): "))

results = []
cutoff_results = []

os.makedirs(CUTOFF_FOLDER, exist_ok=True)

ligand_folders = [
    f for f in os.listdir(DOCKING_FOLDER)
    if os.path.isdir(os.path.join(DOCKING_FOLDER, f))
]

for ligand in sorted(ligand_folders):

    ligand_folder = os.path.join(DOCKING_FOLDER, ligand)
    log_path = os.path.join(ligand_folder, "log.txt")

    if not os.path.exists(log_path):
        continue

    energy = None

    with open(log_path, "r", encoding="utf-8", errors="ignore") as file:
        for line in file:

            line = line.strip()

            if line.startswith("1 "):
                parts = line.split()

                if len(parts) >= 2:
                    try:
                        energy = float(parts[1])
                    except ValueError:
                        pass

                break

    if energy is None:
        continue

    results.append((ligand, energy))

    if energy <= cutoff:

        cutoff_results.append((ligand, energy))

        destination = os.path.join(CUTOFF_FOLDER, ligand)

        if os.path.exists(destination):
            shutil.rmtree(destination)

        shutil.copytree(ligand_folder, destination)

with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as file:

    writer = csv.writer(file)

    writer.writerow([
        "Ligand",
        "BindingEnergy(kcal/mol)"
    ])

    for ligand, energy in results:
        writer.writerow([
            ligand,
            energy
        ])

with open(CUTOFF_CSV, "w", newline="", encoding="utf-8") as file:

    writer = csv.writer(file)

    writer.writerow([
        "Ligand",
        "BindingEnergy(kcal/mol)"
    ])

    for ligand, energy in cutoff_results:
        writer.writerow([
            ligand,
            energy
        ])

print("\n" + "=" * 60)
print("Binding Energy Extraction Completed")
print("=" * 60)

print(f"Total ligand folders found : {len(ligand_folders)}")
print(f"Binding energies extracted  : {len(results)}")
print(f"Cutoff                      : {cutoff} kcal/mol")
print(f"Ligands passing cutoff      : {len(cutoff_results)}")

print("\nFiles created:")
print(f"All binding energies : {OUTPUT_CSV}")
print(f"Cutoff results       : {CUTOFF_CSV}")
print(f"Cutoff folders       : {CUTOFF_FOLDER}")