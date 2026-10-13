import os
import re
import math
import csv
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DOCKING_RESULTS = PROJECT_DIR / "Docking_Results"
OUTPUT_FOLDER = PROJECT_DIR / "Best_Poses"
SUMMARY_FILE = PROJECT_DIR / "best_pose_summary.csv"
INTERACTION_FILE = PROJECT_DIR / "interaction_summary.csv"

SCREENED_FOLDER = PROJECT_DIR / "Screened_top10"
TOP10_FILE = PROJECT_DIR / "top10_screened.csv"

DISTANCE_CUTOFF = 3.0
HYDROGEN_BOND_CUTOFF = 3.0
SALT_BRIDGE_CUTOFF = 4.0
PI_DISTANCE_CUTOFF = 5.5
PI_ANGLE_CUTOFF = 30.0
CATION_PI_DISTANCE_CUTOFF = 6.0
HALOGEN_BOND_CUTOFF = 4.0

HBOND_WEIGHT = 4
SALT_BRIDGE_WEIGHT = 4
PI_WEIGHT = 3
HYDROPHOBIC_WEIGHT = 1

print("=" * 70)
print("BEST POSE SELECTION AND INTERACTION ANALYSIS")
print("=" * 70)

protein_name = input("\nEnter protein PDB file name: ").strip()
PROTEIN_FILE = PROJECT_DIR / protein_name

if not PROTEIN_FILE.exists():
    raise FileNotFoundError(f"\nProtein file not found:\n{PROTEIN_FILE}")

residue_file = input("\nEnter residue file name\n (Press Enter for residues.txt): ").strip()

if residue_file == "":
    residue_file = "residues.txt"

TARGET_RESIDUES_FILE = PROJECT_DIR / residue_file

if not TARGET_RESIDUES_FILE.exists():
    raise FileNotFoundError(f"\nResidue file not found:\n{TARGET_RESIDUES_FILE}")

if not DOCKING_RESULTS.exists():
    raise FileNotFoundError(f"\nDocking_Results folder not found:\n{DOCKING_RESULTS}")

OUTPUT_FOLDER.mkdir(exist_ok=True)

class Atom:
    def __init__(self, atom_name, residue_name, chain, residue_number, x, y, z, element="", atom_type="", charge=0.0):
        self.atom_name = atom_name
        self.residue_name = residue_name
        self.chain = chain
        self.residue_number = residue_number
        self.x = x
        self.y = y
        self.z = z
        self.element = element.upper().strip()
        self.atom_type = atom_type.upper().strip()
        self.charge = charge

def distance(atom1, atom2):
    return math.sqrt(
        (atom1.x - atom2.x) ** 2 +
        (atom1.y - atom2.y) ** 2 +
        (atom1.z - atom2.z) ** 2
    )

def read_target_residues(filename):
    residues = set()
    with open(filename, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            line = line.replace(":", " ")
            fields = line.split()

            if len(fields) >= 2:
                chain = fields[0]

                try:
                    residue = int(fields[1])
                except ValueError:
                    continue

                residues.add((chain, residue))
    return residues

def read_protein(filename):
    atoms = []
    with open(filename, "r", encoding="utf-8", errors="ignore") as file:
        for line in file:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            try:
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain = line[21].strip()
                residue_number = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[76:78].strip()
                atoms.append(Atom(atom_name, residue_name, chain, residue_number, x, y, z, element))
            except (ValueError, IndexError):
                continue
    return atoms

def parse_pose_blocks(filename):
    poses = []
    with open(filename, "r", encoding="utf-8", errors="ignore") as file:
        lines = file.readlines()

    current_pose = []
    pose_number = 0
    inside_model = False
    for line in lines:
        if line.startswith("MODEL"):
            inside_model = True
            current_pose = []
            pose_number += 1
            current_pose.append(line)
            continue

        if inside_model:
            current_pose.append(line)

        if line.startswith("ENDMDL") and inside_model:
            poses.append({"number": pose_number, "lines": current_pose})
            current_pose = []
            inside_model = False

    if len(poses) == 0:
        atom_lines = [
            line
            for line in lines
            if line.startswith(("ATOM", "HETATM"))
        ]

        if atom_lines:
            poses.append({"number": 1, "lines": atom_lines})
    return poses

def extract_binding_energy(lines):
    for line in lines:
        if "REMARK VINA RESULT:" in line:
            match = re.search(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)", line)
            if match:
                return float(match.group(1))
    return None

def get_element_from_atom_type(atom_type, atom_name):
    atom_type = atom_type.upper().strip()
    if atom_type:
        if atom_type.startswith(("CL", "BR")):
            return atom_type[:2]
        if atom_type[0] in "FCINOS":
            return atom_type[0]

    atom_name = re.sub(r"[^A-Za-z]", "", atom_name.upper())
    if atom_name.startswith("CL"):
        return "CL"
    if atom_name.startswith("BR"):
        return "BR"
    if atom_name:
        return atom_name[0]

    return ""

def read_ligand_from_pose(lines):
    atoms = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue

        try:
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            charge = 0.0

            try:
                charge_text = line[70:76].strip()
                if charge_text:
                    charge = float(charge_text)
            except ValueError:
                charge = 0.0

            atom_type = ""
            if len(line) >= 79:
                atom_type = line[77:79].strip()

            element = get_element_from_atom_type(atom_type, atom_name)
            atoms.append(Atom(atom_name, "LIG", "L", 1, x, y, z, element, atom_type, charge))

        except (ValueError, IndexError):
            continue

    return atoms

def get_interacting_residues(protein_atoms, ligand_atoms):
    interacting = set()
    for protein_atom in protein_atoms:
        for ligand_atom in ligand_atoms:
            d = distance (protein_atom, ligand_atom)

            if d <= DISTANCE_CUTOFF:
                interacting.add((protein_atom.chain, protein_atom.residue_number))
                break

    return interacting

def count_required_matches(interacting, required_residues):
    return interacting.intersection(required_residues)

def residue_label(atom):
    return f"{atom.residue_name}{atom.residue_number}({atom.chain})"

def protein_residue_atoms(protein_atoms, target_residues):
    grouped = {}
    for atom in protein_atoms:
        key = (atom.chain, atom.residue_number)
        if key not in target_residues:
            continue

        if key not in grouped:
            grouped[key] = []
        grouped[key].append(atom)

    return grouped

def is_protein_hbond_donor(atom):
    residue = atom.residue_name.upper()
    name = atom.atom_name.upper()
    if atom.element == "N":
        if residue in {"PRO"} and name == "N":
            return False
        return True

    if atom.element == "O":
        if name in {"OH", "OG", "OG1", "OH1"}:
            return True

    if atom.element == "S":
        if name in {"SG", "SD"}:
            return True

    return False

def is_protein_hbond_acceptor(atom):
    residue = atom.residue_name.upper()
    name = atom.atom_name.upper()
    if atom.element == "O":
        if residue == "SER" and name == "OG":
            return True
        if residue == "THR" and name == "OG1":
            return True
        if residue == "TYR" and name == "OH":
            return True
        return True

    if atom.element == "N":
        if residue == "LYS":
            return False
        if residue == "ARG":
            return False
        if residue == "PRO" and name == "N":
            return False
        if residue == "HIS" and name in {"ND1", "NE2"}:
            return True
        if residue in {"ASN", "GLN"} and name in {"OD1", "OE1"}:
            return False
        return True

    if atom.element == "S":
        if residue in {"CYS", "MET"}:
            return True

    return False

def is_ligand_hbond_donor(atom):
    if atom.element == "N":
        return True
    if atom.element == "O":
        return True
    if atom.element == "S":
        return True
    return False


def is_ligand_hbond_acceptor(atom):
    if atom.element in {"O", "S"}:
        return True
    if atom.element == "N":
        atom_type = atom.atom_type.upper()
        if atom_type in {"NA", "N", "N3", "N2", "N1"}:
            return True
        if atom.charge > 0.5:
            return False

        return True
    return False

def is_positive_protein_atom(atom):
    residue = atom.residue_name.upper()
    name = atom.atom_name.upper()
    if residue == "LYS" and name == "NZ":
        return True
    if residue == "ARG" and name in {"NE", "NH1", "NH2"}:
        return True
    if residue == "HIS" and name in {"ND1", "NE2"}:
        return True

    return False

def is_negative_protein_atom(atom):
    residue = atom.residue_name.upper()
    name = atom.atom_name.upper()
    if residue == "ASP" and name in {"OD1", "OD2"}:
        return True
    if residue == "GLU" and name in {"OE1", "OE2"}:
        return True

    return False

def is_positive_ligand_atom(atom):
    if atom.charge > 0.5:
        return True
    if atom.element == "N":
        if atom.atom_type.upper() in {"N3", "N4", "NA"}:
            return True

    return False

def is_negative_ligand_atom(atom):
    if atom.charge < -0.5:
        return True
    if atom.element == "O":
        if atom.atom_type.upper() in {"OA", "O2"}:
            return True

    return False

def is_hydrophobic_residue(residue_name):
    return residue_name.upper() in {"ALA", "VAL", "LEU", "ILE", "MET", "PRO", "PHE", "TRP", "TYR"}

def is_ligand_hydrophobic_atom(atom):
    return atom.element in {"C", "S"}

def get_aromatic_atoms(atoms):
    aromatic = []
    for atom in atoms:
        atom_type = atom.atom_type.upper()
        if atom_type in {"A", "C", "N"}:
            if atom.element in {"C", "N"}:
                aromatic.append(atom)

    return aromatic

def centroid(atoms):
    if not atoms:
        return None
    x = sum(atom.x for atom in atoms) / len(atoms)
    y = sum(atom.y for atom in atoms) / len(atoms)
    z = sum(atom.z for atom in atoms) / len(atoms)
    return (x, y, z)

def centroid_distance(c1, c2):
    if c1 is None or c2 is None:
        return float("inf")
    return math.sqrt(
        (c1[0] - c2[0]) ** 2 +
        (c1[1] - c2[1]) ** 2 +
        (c1[2] - c2[2]) ** 2
    )

def ring_normal(atoms):
    if len(atoms) < 3:
        return None
    p1 = atoms[0]
    p2 = atoms[1]
    p3 = atoms[2]

    v1 = (p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
    v2 = (p3.x - p1.x, p3.y - p1.y, p3.z - p1.z)
    nx = v1[1] * v2[2] - v1[2] * v2[1]
    ny = v1[2] * v2[0] - v1[0] * v2[2]
    nz = v1[0] * v2[1] - v1[1] * v2[0]

    norm = math.sqrt(nx * nx + ny * ny + nz * nz)

    if norm == 0:
        return None

    return (nx / norm, ny / norm, nz / norm)

def vector_between(c1, c2):
    return (c2[0] - c1[0], c2[1] - c1[1], c2[2] - c1[2])

def vector_angle(v1, v2):
    dot = (
        v1[0] * v2[0] +
        v1[1] * v2[1] +
        v1[2] * v2[2]
    )

    n1 = math.sqrt(
        v1[0] ** 2 +
        v1[1] ** 2 +
        v1[2] ** 2
    )

    n2 = math.sqrt(
        v2[0] ** 2 +
        v2[1] ** 2 +
        v2[2] ** 2
    )

    if n1 == 0 or n2 == 0:
        return 180.0

    value = dot / (n1 * n2)
    value = max(-1.0, min(1.0, value))
    return math.degrees(math.acos(value))

def get_aromatic_groups(atoms):
    groups = []
    aromatic_atoms = get_aromatic_atoms(atoms)
    if len(aromatic_atoms) >= 5:
        groups.append(aromatic_atoms)
    return groups

def format_residue_with_distance(atom, distance_value):
    return f"{residue_label(atom)} [{distance_value:.2f} Å]"

def add_interaction(dictionary, key, atom, distance_value):
    value = format_residue_with_distance(atom, distance_value)
    if key not in dictionary:
        dictionary[key] = set()
    dictionary[key].add(value)

def analyze_interactions(protein_atoms, ligand_atoms, target_residues):
    interaction_data = {"Hydrogen bonds": set(), "Hydrophobic interactions": set(), "Salt bridges": set(), "π–π interactions": set(), "Cation–π interactions": set(), "Halogen bonds": set()}
    grouped_protein = protein_residue_atoms(protein_atoms, target_residues)
    for residue_key, residue_atoms in grouped_protein.items():
        if not residue_atoms:
            continue
        residue_name = residue_atoms[0].residue_name

        for protein_atom in residue_atoms:
            for ligand_atom in ligand_atoms:
                d = distance(protein_atom, ligand_atom)
                if d <= HYDROGEN_BOND_CUTOFF:
                    if (is_protein_hbond_donor(protein_atom) and is_ligand_hbond_acceptor(ligand_atom)):
                        add_interaction(interaction_data, "Hydrogen bonds", protein_atom, d)

                    elif (is_ligand_hbond_donor(ligand_atom) and is_protein_hbond_acceptor(protein_atom)):
                        add_interaction(interaction_data, "Hydrogen bonds", protein_atom, d)

                if d <= SALT_BRIDGE_CUTOFF:
                    if (is_positive_protein_atom(protein_atom) and is_negative_ligand_atom(ligand_atom)):
                        add_interaction(interaction_data, "Salt bridges", protein_atom, d)

                    elif (is_negative_protein_atom(protein_atom) and is_positive_ligand_atom(ligand_atom)):
                        add_interaction(interaction_data, "Salt bridges", protein_atom, d)

                if d <= DISTANCE_CUTOFF:
                    if (is_hydrophobic_residue(residue_name) and is_ligand_hydrophobic_atom(ligand_atom)):
                        add_interaction(interaction_data, "Hydrophobic interactions", protein_atom, d)

                if d <= HALOGEN_BOND_CUTOFF:
                    if (ligand_atom.element in {"F", "CL", "BR", "I"} and protein_atom.element in {"O", "N", "S"}):
                        add_interaction(interaction_data, "Halogen bonds", protein_atom, d)
    protein_aromatic_groups = []

    for residue_key, residue_atoms in grouped_protein.items():
        residue_name = residue_atoms[0].residue_name.upper()

        if residue_name in {"PHE", "TYR", "TRP", "HIS"}:
            group = get_aromatic_groups(residue_atoms)

            if group:
                protein_aromatic_groups.append((residue_key, group[0]))

    ligand_aromatic_groups = get_aromatic_groups(ligand_atoms)

    for residue_key, protein_ring in protein_aromatic_groups:
        p_centroid = centroid(protein_ring)

        for ligand_ring in ligand_aromatic_groups:
            l_centroid = centroid(ligand_ring)
            d = centroid_distance(p_centroid, l_centroid)
            if d <= PI_DISTANCE_CUTOFF:
                p_normal = ring_normal(protein_ring)
                l_normal = ring_normal(ligand_ring)
                valid_geometry = True

                if p_normal and l_normal:
                    angle = vector_angle(p_normal, l_normal)
                    parallel_angle = min(angle, abs(180.0 - angle))
                    valid_geometry = (parallel_angle <= PI_ANGLE_CUTOFF)

                if valid_geometry:
                    residue_atom = protein_ring[0]
                    add_interaction(interaction_data, "π–π interactions", residue_atom, d)
    protein_aromatic_groups_for_cation = []

    for residue_key, residue_atoms in grouped_protein.items():
        residue_name = residue_atoms[0].residue_name.upper()
        if residue_name in {"PHE", "TYR", "TRP", "HIS"}:
            group = get_aromatic_groups(residue_atoms)
            if group:
                protein_aromatic_groups_for_cation.append((residue_key, group[0]))

    for residue_key, protein_ring in protein_aromatic_groups_for_cation:
        p_centroid = centroid(protein_ring)
        for ligand_atom in ligand_atoms:
            if not is_positive_ligand_atom(ligand_atom):
                continue
            d = centroid_distance(p_centroid, (ligand_atom.x, ligand_atom.y, ligand_atom.z))
            if d <= CATION_PI_DISTANCE_CUTOFF:
                residue_atom = protein_ring[0]
                add_interaction(interaction_data, "Cation–π interactions", residue_atom, d)

    for residue_key, residue_atoms in grouped_protein.items():
        aromatic = get_aromatic_groups(residue_atoms)
        if not aromatic:
            continue

        ring = aromatic[0]
        ring_centroid = centroid(ring)
        for ligand_atom in ligand_atoms:
            if not is_positive_protein_atom(ligand_atom):
                continue
            d = centroid_distance(ring_centroid, (ligand_atom.x, ligand_atom.y, ligand_atom.z))
            if d <= CATION_PI_DISTANCE_CUTOFF:
                residue_atom = ring[0]
                add_interaction(interaction_data, "Cation–π interactions", residue_atom, d)
    return interaction_data

def clean_interaction_values(values):
    if not values:
        return ""
    return "; ".join(sorted(values))

def select_best_pose(ligand_folder, protein_atoms, required_residues):
    out_file = ligand_folder / "out.pdbqt"
    if not out_file.exists():
        return None, "No out.pdbqt", None
    if out_file.stat().st_size == 0:
        return None, "Empty out.pdbqt", None
    poses = parse_pose_blocks(out_file)
    if not poses:
        return None, "No poses found", None
    pose_data = []

    for pose in poses:
        ligand_atoms = read_ligand_from_pose(pose["lines"])
        if not ligand_atoms:
            continue

        interacting = get_interacting_residues(protein_atoms, ligand_atoms)
        matches = count_required_matches(interacting, required_residues)
        energy = extract_binding_energy(pose["lines"])
        interactions = analyze_interactions(protein_atoms, ligand_atoms, required_residues)
        pose_data.append({"number": pose["number"], "lines": pose["lines"], "energy": energy, "interacting": interacting, "matches": matches, "match_count": len(matches), "interactions": interactions})

    if not pose_data:
        return None, "No valid poses", None

    matching_poses = [
        pose
        for pose in pose_data
        if pose["match_count"] >= 1
    ]

    if not matching_poses:
        return (None, "No pose interacts with target residues", len(poses))

    matching_poses.sort(
        key=lambda x: (
            -x["match_count"],
            x["energy"]
            if x["energy"] is not None
            else float("inf")
        )
    )

    best = matching_poses[0]

    return (best, "Target residue interaction found", len(poses))

def convert_pose_to_pdb(lines):
    pdb_lines = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            pdb_line = line.rstrip()
            if len(pdb_line) >= 66:
                pdb_line = pdb_line[:66]
            pdb_lines.append(pdb_line + "\n")
    pdb_lines.append("END\n")
    return pdb_lines

def process_all_ligands():
    protein_atoms = read_protein(PROTEIN_FILE)
    required_residues = read_target_residues(TARGET_RESIDUES_FILE)
    if not protein_atoms:
        raise RuntimeError("No protein atoms were found.")
    if not required_residues:
        raise RuntimeError("No target residues were found in residues.txt.")

    print(f"\nProtein atoms read: {len(protein_atoms)}")
    print(f"Target residues read: {len(required_residues)}")
    print("\nTarget residues:")
    for chain, residue in sorted(required_residues):
        print(f"  {chain}:{residue}")

    folders = sorted(
        [
            folder
            for folder in DOCKING_RESULTS.iterdir()
            if folder.is_dir()
        ]
    )

    summary = []
    interaction_summary = []

    print("\n")
    print("=" * 70)
    print(f"Docking folders found: {len(folders)}")
    print("=" * 70)

    for number, folder in enumerate(folders, start=1):
        print(f"\n[{number}/{len(folders)}] Processing: {folder.name}")
        best_pose, reason, total_poses = (select_best_pose (folder, protein_atoms, required_residues))
        if best_pose is None:
            print(f"Skipped: {reason}")
            summary.append([folder.name, "", "", "", "", reason, ""])
            continue

        pose_number = best_pose["number"]
        energy = best_pose["energy"]
        match_count = best_pose["match_count"]
        output_name = (f"{folder.name}_pose{pose_number}.pdb")
        destination = (OUTPUT_FOLDER / output_name)
        pdb_lines = convert_pose_to_pdb(best_pose["lines"])

        with open(destination, "w", encoding="utf-8") as file:
            file.writelines(pdb_lines)

        interacting_residues = "; ".join(
            [
                f"{chain}:{residue}"
                for chain, residue
                in sorted(best_pose["interacting"])
            ]
        )

        matched_residues = "; ".join(
            [
                f"{chain}:{residue}"
                for chain, residue
                in sorted(best_pose["matches"])
            ]
        )

        interaction_data = best_pose["interactions"]
        interaction_row = [folder.name, pose_number, energy if energy is not None else "", matched_residues, clean_interaction_values (interaction_data ["Hydrogen bonds"]),
            clean_interaction_values(interaction_data ["Hydrophobic interactions"]), clean_interaction_values (interaction_data["Salt bridges"]),
            clean_interaction_values(interaction_data ["π–π interactions"]), clean_interaction_values (interaction_data ["Cation–π interactions"]),
            clean_interaction_values(interaction_data ["Halogen bonds"]), interacting_residues, output_name]

        interaction_summary.append(interaction_row)
        print(f"Total poses: {total_poses}")
        print(f"Selected pose: {pose_number}")
        print(
            f"Binding energy: "
            f"{energy}"
            if energy is not None
            else
            "Binding energy: Not found"
        )

        print(f"Target residue matches: {match_count}")
        print(f"Matched residues: {matched_residues}")
        for interaction_type in ["Hydrogen bonds", "Hydrophobic interactions", "Salt bridges", "π–π interactions", "Cation–π interactions", "Halogen bonds"]:
            values = clean_interaction_values(interaction_data [interaction_type])
            print(f"{interaction_type}: {values if values else 'None'}")
        print(f"Saved: {destination.name}")
        summary.append([folder.name, pose_number, energy if energy is not None else "", matched_residues, interacting_residues, reason, destination.name])

    with open(SUMMARY_FILE, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Ligand", "Selected_Pose", "Binding_Energy_kcal_mol", "Matched_Required_Residues", "All_Interacting_Residues", "Selection_Reason", "Output_File"])
        writer.writerows(summary)

    with open(INTERACTION_FILE, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Ligand", "Selected_Pose", "Binding_Energy_kcal_mol", "Matched_Required_Residues", "Hydrogen bonds", "Hydrophobic interactions",
        "Salt bridges", "π–π interactions", "Cation–π interactions", "Halogen bonds", "All_Interacting_Residues", "Output_File"])
        writer.writerows(interaction_summary)

    print("\n")
    print("=" * 70)
    print("BEST POSE AND INTERACTION ANALYSIS COMPLETED")
    print("=" * 70)
    print(f"Best poses saved in : {OUTPUT_FOLDER}")
    print(f"Pose summary saved in : {SUMMARY_FILE}")
    print(f"Interaction summary saved in : {INTERACTION_FILE}")
    print(f"Number of selected molecules: {len(interaction_summary)}")
    print("=" * 70)

def count_interactions(value):
    if not value or not str(value).strip():
        return 0
    return len([x for x in str(value).split(";") if x.strip()])

def screen_top_10_percent():
    if not INTERACTION_FILE.exists():
        print("\nInteraction summary file not found.")
        return

    with open(INTERACTION_FILE, "r", newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)

    if not rows:
        print("\nNo molecules found in interaction summary.")
        return

    processed_rows = []

    for row in rows:
        try:
            energy = float(row["Binding_Energy_kcal_mol"])
        except (ValueError, TypeError):
            continue

        hbond_count = count_interactions(row["Hydrogen bonds"])
        salt_bridge_count = count_interactions(row["Salt bridges"])
        pi_count = count_interactions(row["π–π interactions"])
        hydrophobic_count = count_interactions(row["Hydrophobic interactions"])

        interaction_score = (HBOND_WEIGHT * hbond_count + SALT_BRIDGE_WEIGHT * salt_bridge_count + PI_WEIGHT * pi_count + HYDROPHOBIC_WEIGHT * hydrophobic_count)
        row["Hydrogen_Bond_Count"] = hbond_count
        row["Salt_Bridge_Count"] = salt_bridge_count
        row["Pi_Interaction_Count"] = pi_count
        row["Hydrophobic_Interaction_Count"] = hydrophobic_count
        row["Interaction_Score"] = interaction_score
        processed_rows.append(row)

    if not processed_rows:
        print("\nNo valid molecules available for screening.")
        return

    processed_rows.sort(key=lambda x: (float(x["Binding_Energy_kcal_mol"]), -float(x["Interaction_Score"]), -int(x["Hydrogen_Bond_Count"]), -int(x["Salt_Bridge_Count"]), -int(x["Pi_Interaction_Count"]),
            -int(x["Hydrophobic_Interaction_Count"])))

    top_n = max(1, math.ceil(len(processed_rows) * 0.10))
    top10_rows = processed_rows[:top_n]

    SCREENED_FOLDER.mkdir(exist_ok=True)
    for file in SCREENED_FOLDER.iterdir():
        if file.is_file():
            file.unlink()

    output_columns = ["Ligand", "Selected_Pose", "Binding_Energy_kcal_mol", "Matched_Required_Residues", "Hydrogen bonds", "Hydrogen_Bond_Count", "Salt bridges",
        "Salt_Bridge_Count", "π–π interactions", "Pi_Interaction_Count", "Hydrophobic interactions", "Hydrophobic_Interaction_Count", "Cation–π interactions",
        "Halogen bonds", "Interaction_Score", "All_Interacting_Residues", "Output_File"]

    with open(TOP10_FILE, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=output_columns, extrasaction="ignore")
        writer.writeheader()
        for row in top10_rows:
            writer.writerow(row)
    copied = 0

    for row in top10_rows:
        output_file = row["Output_File"]

        if not output_file:
            continue

        source = OUTPUT_FOLDER / output_file
        destination = SCREENED_FOLDER / output_file
        if source.exists():
            with open(source, "rb") as source_file:
                with open(destination, "wb") as destination_file:
                    destination_file.write(source_file.read())

            copied += 1

    print("\n")
    print("=" * 70)
    print("TOP 10% SCREENING COMPLETED")
    print("=" * 70)
    print(f"Total valid molecules : {len(processed_rows)}")
    print(f"Top 10% selected      : {top_n}")
    print(f"Molecules copied      : {copied}")
    print(f"Screened folder       : {SCREENED_FOLDER}")
    print(f"Top 10% CSV            : {TOP10_FILE}")
    print("=" * 70)
    print("\nTop selected molecules:")

    for rank, row in enumerate(top10_rows, start=1):
        print(
            f"{rank:>3}. "
            f"{row['Ligand']} | "
            f"Energy: {row['Binding_Energy_kcal_mol']} | "
            f"HB: {row['Hydrogen_Bond_Count']} | "
            f"SB: {row['Salt_Bridge_Count']} | "
            f"π: {row['Pi_Interaction_Count']} | "
            f"Hydrophobic: {row['Hydrophobic_Interaction_Count']} | "
            f"Score: {row['Interaction_Score']}"
        )

def main():
    process_all_ligands()
    screen_top_10_percent()

if __name__ == "__main__":
    main()