import sys
import csv
import math
from pathlib import Path
from collections import defaultdict

AMINO_ACIDS = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
EXCLUDED_HETERO = {"HOH", "WAT", "DOD", "SOL", "NA", "CL", "MG", "MN", "ZN", "CA", "FE", "K", "CU", "CO", "NI", "CD"}

GRID_PADDING = 5.0
INTERACTION_CUTOFF = 4.0

def parse_pdb(filename):
    atoms = []
    with open(filename, "r", encoding="utf-8", errors="ignore") as file:
        for line in file:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            try:
                record = line[0:6].strip()
                atom_name = line[12:16].strip()
                residue = line[17:20].strip().upper()
                chain = line[21].strip()
                residue_number = int(line[22:26].strip())

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                element = line[76:78].strip()
                if not element:
                    element = atom_name[0]
                atoms.append({"record": record, "atom":atom_name, "residue":residue, "chain":chain, "resnum":residue_number, "x":x, "y":y, "z":z, "element":element.upper()})

            except (ValueError, IndexError):
                continue

    return atoms


def identify_ligands(atoms):
    ligand_groups = defaultdict(list)
    for atom in atoms:
        if atom["record"] != "HETATM":
            continue
        residue = atom["residue"]
        if residue in EXCLUDED_HETERO:
            continue
        if residue in AMINO_ACIDS:
            continue
        key = (residue, atom["chain"], atom["resnum"])
        ligand_groups[key].append(atom)
    candidates = []

    for key, group in ligand_groups.items():
        heavy_atoms = [
            atom for atom in group
            if atom["element"] != "H"
        ]

        if len(heavy_atoms) >= 3:
            candidates.append((key, group))

    candidates.sort(key=lambda x: len(x[1]), reverse=True)
    return candidates

def select_ligand(atoms):
    candidates = identify_ligands(atoms)
    if not candidates:
        raise RuntimeError("No organic ligand could be identified from the PDB file.")
    print("\nLigand candidates found:")
    for number, (key, group) in enumerate(candidates, start=1):
        residue, chain, resnum = key
        heavy_atoms = sum(
            1
            for atom in group
            if atom["element"] != "H"
        )

        print(
            f"{number}. "
            f"{residue} "
            f"chain {chain or '-'} "
            f"residue {resnum} "
            f"({heavy_atoms} heavy atoms)"
        )

    if len(candidates) == 1:
        selected = candidates[0]
        print(
            f"\nAutomatically selected ligand: "
            f"{selected[0][0]} "
            f"chain {selected[0][1] or '-'} "
            f"residue {selected[0][2]}"
        )
        return selected[1]

    while True:
        try:
            choice = input("\nEnter the number of the ligand to analyze: ").strip()
            number = int(choice)
            if 1 <= number <= len(candidates):
                selected = candidates[number - 1]
                print(
                    f"\nSelected ligand: "
                    f"{selected[0][0]} "
                    f"chain {selected[0][1] or '-'} "
                    f"residue {selected[0][2]}"
                )
                return selected[1]
        except ValueError:
            pass
        print("Invalid selection. Try again.")

def get_protein_atoms(atoms):
    return [
        atom
        for atom in atoms
        if (
            atom["record"] == "ATOM" and
            atom["residue"] in AMINO_ACIDS
        )
    ]

def calculate_distance(atom1, atom2):
    return math.sqrt(
        (atom1["x"] - atom2["x"]) ** 2 +
        (atom1["y"] - atom2["y"]) ** 2 +
        (atom1["z"] - atom2["z"]) ** 2
    )

def calculate_gridbox(ligand_atoms):
    x_values = [
        atom["x"]
        for atom in ligand_atoms
    ]

    y_values = [
        atom["y"]
        for atom in ligand_atoms
    ]

    z_values = [
        atom["z"]
        for atom in ligand_atoms
    ]

    min_x = min(x_values)
    max_x = max(x_values)

    min_y = min(y_values)
    max_y = max(y_values)

    min_z = min(z_values)
    max_z = max(z_values)

    center_x = (min_x + max_x) / 2.0
    center_y = (min_y + max_y) / 2.0
    center_z = (min_z + max_z) / 2.0

    size_x = (max_x - min_x) + (2 * GRID_PADDING)
    size_y = (max_y - min_y) + (2 * GRID_PADDING)
    size_z = (max_z - min_z) + (2 * GRID_PADDING)

    return {"center_x": center_x, "center_y": center_y, "center_z": center_z, "size_x": size_x, "size_y": size_y, "size_z": size_z, "min_x": min_x,
        "max_x": max_x, "min_y": min_y, "max_y": max_y, "min_z": min_z, "max_z": max_z}

def find_interacting_residues(ligand_atoms, protein_atoms):
    residue_contacts = defaultdict(lambda: {"minimum_distance": float("inf"), "ligand_atoms": set(), "protein_atoms": set()})
    for ligand_atom in ligand_atoms:
        for protein_atom in protein_atoms:
            d = calculate_distance(ligand_atom, protein_atom)
            if d <= INTERACTION_CUTOFF:
                key = (protein_atom["chain"], protein_atom["resnum"], protein_atom["residue"])

                if d < residue_contacts[key]["minimum_distance"]:
                    residue_contacts[key]["minimum_distance"] = d

                residue_contacts[key]["ligand_atoms"].add(ligand_atom["atom"])
                residue_contacts[key]["protein_atoms"].add(protein_atom["atom"])

    results = []
    for key, data in residue_contacts.items():
        chain, resnum, residue = key
        results.append({"chain": chain, "resnum": resnum, "residue": residue, "minimum_distance": data["minimum_distance"],
            "ligand_atoms": sorted(data["ligand_atoms"]), "protein_atoms": sorted(data["protein_atoms"])})
    results.sort(key=lambda x: x["minimum_distance"])
    return results

def write_gridbox_file(filename, grid, ligand_atoms):
    with open(filename, "w", encoding="utf-8") as file:
        file.write("AutoDock Vina Grid Box\n")
        file.write(f"center_x = {grid['center_x']:.3f}\n")
        file.write(f"center_y = {grid['center_y']:.3f}\n")
        file.write(f"center_z = {grid['center_z']:.3f}\n")
        file.write("\n")
        file.write(f"size_x = {grid['size_x']:.3f}\n")
        file.write(f"size_y = {grid['size_y']:.3f}\n")
        file.write(f"size_z = {grid['size_z']:.3f}\n")
        file.write("\n")
        file.write(f"Padding = {GRID_PADDING:.1f} Angstrom\n")
        file.write("\n")
        file.write("Ligand coordinate range\n")
        file.write("\n")
        file.write(f"X: {grid['min_x']:.3f} to {grid['max_x']:.3f}\n")
        file.write(f"Y: {grid['min_y']:.3f} to {grid['max_y']:.3f}\n")
        file.write(f"Z: {grid['min_z']:.3f} to {grid['max_z']:.3f}\n")
        file.write("\n")
        file.write(f"Ligand atoms: {len(ligand_atoms)}\n")

def write_vina_config(filename, grid):
    with open(filename, "w", encoding="utf-8") as file:
        file.write(f"center_x = {grid['center_x']:.3f}\n")
        file.write(f"center_y = {grid['center_y']:.3f}\n")
        file.write(f"center_z = {grid['center_z']:.3f}\n")
        file.write("\n")
        file.write(f"size_x = {grid['size_x']:.3f}\n")
        file.write(f"size_y = {grid['size_y']:.3f}\n")
        file.write(f"size_z = {grid['size_z']:.3f}\n")

def write_interactions(filename, interactions):
    with open(filename, "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["Chain", "Residue Number", "Residue", "Minimum Distance (Angstrom)", "Protein Atoms", "Ligand Atoms"])
        for item in interactions:
            writer.writerow([item["chain"], item["resnum"], item["residue"], f"{item['minimum_distance']:.3f}", "; ".join( item["protein_atoms"]), "; ".join(item["ligand_atoms"])])

def write_residues_file(filename, interactions):
    with open(filename, "w", encoding="utf-8") as file:
        for item in interactions:
            if item["chain"]:
                file.write(f"{item['chain']}: {item['resnum']}\n")
            else:
                file.write(f"{item['resnum']}\n")

def main():
    if len(sys.argv) != 2:
        print("Usage:")
        print("python gridbox_interactions.py protein_ligand.pdb")
        return

    pdb_file = Path(sys.argv[1]).resolve()
    if not pdb_file.exists():
        print(f"File not found: {pdb_file}")
        return

    output_folder = (
        pdb_file.parent /
        f"{pdb_file.stem}_analysis"
    )

    output_folder.mkdir(exist_ok=True)

    gridbox_file = (
        output_folder /
        "gridbox.txt"
    )

    vina_config_file = (
        output_folder /
        "vina_config.txt"
    )

    interactions_file = (
        output_folder /
        "interacting_residues.csv"
    )

    residues_file = (
        output_folder /
        "residues.txt"
    )

    print("=" * 65)
    print("PROTEIN-LIGAND GRID BOX AND INTERACTION ANALYSIS")
    print("=" * 65)

    print(f"\nInput file: {pdb_file.name}")
    atoms = parse_pdb(pdb_file)
    print(f"Total atoms: {len(atoms)}")
    ligand_atoms = select_ligand(atoms)
    protein_atoms = get_protein_atoms(atoms)
    print(f"Protein atoms: {len(protein_atoms)}")
    print(f"Ligand atoms: {len(ligand_atoms)}")
    grid = calculate_gridbox(ligand_atoms)
    print("\nGRID BOX")
    print("-" * 30)
    print(f"Center X = {grid['center_x']:.3f}")
    print(f"Center Y = {grid['center_y']:.3f}")
    print(f"Center Z = {grid['center_z']:.3f}")
    print(f"Size X   = {grid['size_x']:.3f}")
    print(f"Size Y   = {grid['size_y']:.3f}")
    print(f"Size Z   = {grid['size_z']:.3f}")

    interactions = find_interacting_residues(ligand_atoms, protein_atoms)
    print("\nINTERACTING RESIDUES")
    print("-" * 50)
    for item in interactions:
        chain = (
            item["chain"]
            if item["chain"]
            else "-"
        )

        print(
            f"{chain}:"
            f"{item['residue']}"
            f"{item['resnum']}   "
            f"{item['minimum_distance']:.3f} Å"
        )

    write_gridbox_file(gridbox_file, grid, ligand_atoms)
    write_vina_config(vina_config_file, grid)
    write_interactions(interactions_file, interactions)
    write_residues_file(residues_file, interactions)

    print("\n")
    print("=" * 65)
    print("ANALYSIS COMPLETED")
    print("=" * 65)
    print(f"\nGrid box information: \n{gridbox_file}")
    print(f"\nVina configuration: \n{vina_config_file}")
    print(f"\nInteracting residues: \n{interactions_file}")
    print(f"\nResidue list: \n{residues_file}")

if __name__ == "__main__":
    main()