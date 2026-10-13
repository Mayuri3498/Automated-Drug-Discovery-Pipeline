import argparse
import os
import re
import shutil
import rdkit
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Image as RLImage, Paragraph, Spacer, PageBreak

NUMERIC_RULES = {
    "MW_<=500": ("MW", "<=", 500),
    "LogP_<=5": ("LogP", "<=", 5),
    "Rotatable_Bonds_<=10": ("Rotatable Bonds", "<=", 10),
    "Acceptors_<=10": ("Acceptors", "<=", 10),
    "Donors_<=5": ("Donors", "<=", 5),
    "Caco2_>=-6": ("Caco-2 (logPaap)", ">=", -6),
    "LogS_>=-5": ("Log S", ">=", -5),
}

CATEGORICAL_RULES = {
    "Human_Intestinal_Absorption": ("Human Intestinal Absorption", "Absorbed"),
    "Human_Oral_Bioavailability_50": ("Human Oral Bioavailability 50%", "Bioavailable"),
    "P_Glycoprotein_Inhibitor": ("P-Glycoprotein Inhibitor", "Non-Inhibitor"),
    "AMES_Mutagenesis": ("AMES Mutagenesis", "Safe"),
    "Carcinogenesis": ("Carcinogenesis", "Safe"),
    "hERG_Blockers": ("hERG Blockers", "Safe"),
    "Liver_Injury_I_DILI": ("Liver Injury I (DILI)", "Safe"),
    "Respiratory_Disease": ("Respiratory Disease", "Safe"),
    "Skin_Sensitisation": ("Skin Sensitisation", "Safe"),
}

PREFIXES = ("Absorption ", "Distribution ", "Metabolism ", "Excretion ", "Toxicity ", "General Properties ",)

def normalize_header(value):
    value = str(value).strip()
    for prefix in PREFIXES:
        if value.startswith(prefix):
            value = value[len(prefix):]
    return value

def normalize_dataframe_columns(df):
    df = df.copy()
    df.rename(columns={c: normalize_header(c) for c in df.columns}, inplace=True)
    return df

def numeric_values(series):
    extracted = (
        series.astype(str)
        .str.replace(",", "", regex=False)
        .str.extract(
            r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
            expand=False,
        )
    )
    return pd.to_numeric(extracted, errors="coerce")

def categorical_values(series):
    return (
        series.astype(str)
        .str.strip()
        .str.replace(r"\s+\(.*\)$", "", regex=True)
        .str.strip()
    )

def find_column(df, requested):
    if requested in df.columns:
        return requested
    target = normalize_header(requested).lower()
    for column in df.columns:
        if normalize_header(column).lower() == target:
            return column
    return None

def clean_identifier(value):
    value = str(value).strip()
    if not value or value.lower() in {"nan", "none"}:
        return ""
    value = os.path.basename(value)
    value = re.sub(
        r"\.(sdf|sd|mol|mol2|pdb|pdbqt|smi|smiles|csv|txt)$",
        "",
        value,
        flags=re.I,
    )
    value = re.sub(r"[^A-Za-z0-9]+", "", value).lower()
    return value

def row_identifiers(row):
    columns = ["Molecule_Label", "Molecule", "Ligand", "Ligand_Name", "Name", "Filename", "File_Name", "Canonical_SMILES", "SMILES", "Molecule_Number",]
    values = []
    for column in columns:
        if column in row.index:
            value = clean_identifier(row[column])
            if value:
                values.append(value)
    return list(dict.fromkeys(values))

def find_pose_file(best_poses_folder, identifiers):
    exact_matches = []
    partial_matches = []

    for root, _, files in os.walk(best_poses_folder):
        for filename in files:
            full_path = os.path.join(root, filename)
            stem = clean_identifier(filename)
            if not stem:
                continue

            if stem in identifiers:
                exact_matches.append(full_path)
                continue

            if any(identifier in stem or stem in identifier
                   for identifier in identifiers if len(identifier) >= 5):
                partial_matches.append(full_path)

    if exact_matches:
        return sorted(exact_matches)[0]

    if partial_matches:
        return sorted(partial_matches)[0]

    return None

def mol_from_structure(path):
    extension = Path(path).suffix.lower()
    if extension == ".sdf":
        supplier = Chem.SDMolSupplier(path, removeHs=False, sanitize=False,)
        mol = next((m for m in supplier if m is not None), None)

    elif extension == ".mol":
        mol = Chem.MolFromMolFile(path, removeHs=False, sanitize=False)

    elif extension == ".mol2":
        mol = Chem.MolFromMol2File(path, removeHs=False, sanitize=False)

    else:
        mol = None

    if mol is not None:
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass

    return mol

def smiles_from_row(row):
    for column in ("Canonical_SMILES", "Canonical SMILES", "SMILES", "Smiles", "smiles",):
        if column in row.index:
            value = str(row[column]).strip()
            if value and value.lower() not in {"nan", "none"}:
                return value
    return ""

def molecule_label(row, fallback):
    for column in ("Molecule_Label", "Molecule", "Ligand", "Ligand_Name", "Name", "Filename", "File_Name", "Molecule_Number",):
        if column in row.index:
            value = str(row[column]).strip()
            if value and value.lower() not in {"nan", "none"}:
                return value
    return f"Molecule_{fallback}"

def safe_filename(value):
    value = re.sub(r'[<>:"/\\|?*]', "_", str(value))
    value = re.sub(r"\s+", "_", value)
    return value[:150] or "Molecule"

def generate_2d_image(mol, output_path):
    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        pass

    image = Draw.MolToImage(mol, size=(1200, 900), kekulize=True,)
    image.save(output_path)

def create_2d_pdf(image_folder, output_pdf, title="Strict-Pass Molecules - 2D Structures"):
    image_folder = Path(image_folder)
    image_files = sorted(
        [p for p in image_folder.glob("*.png") if p.is_file()],
        key=lambda p: p.name,
    )

    if not image_files:
        raise FileNotFoundError(f"No 2D PNG images found in {image_folder}")

    page_width, page_height = landscape(A4)
    doc = SimpleDocTemplate(str(output_pdf), pagesize=landscape(A4), rightMargin=10 * mm, leftMargin=10 * mm, topMargin=10 * mm, bottomMargin=10 * mm,)

    styles = getSampleStyleSheet()
    title_style = ParagraphStyle("MoleculeTitle", parent=styles["Title"], fontSize=14, leading=16, alignment=1, spaceAfter=5 * mm)
    label_style = ParagraphStyle("MoleculeLabel", parent=styles["Normal"], fontSize=7.5, leading=9, alignment=1, wordWrap="CJK",)
    usable_width = page_width - 20 * mm
    col_width = usable_width / 3.0
    usable_height = page_height - 30 * mm
    image_height = 50 * mm
    image_width = col_width - 8 * mm

    story = [Paragraph(title, title_style)]
    cells_per_page = 9
    for start in range(0, len(image_files), cells_per_page):
        batch = image_files[start:start + cells_per_page]
        rows = []
        current_row = []

        for img_path in batch:
            img = RLImage(str(img_path), width=image_width, height=image_height)
            label = img_path.stem.replace("_", " ")
            cell = [img, Spacer(1, 1.5 * mm), Paragraph(label, label_style)]
            current_row.append(cell)

            if len(current_row) == 3:
                rows.append(current_row)
                current_row = []

        if current_row:
            while len(current_row) < 3:
                current_row.append("")
            rows.append(current_row)

        table = Table(rows, colWidths=[col_width] * 3, rowHeights=[58 * mm] * len(rows))
        table.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("ALIGN", (0, 0), (-1, -1), "CENTER"), ("LEFTPADDING", (0, 0), (-1, -1), 2 * mm),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2 * mm), ("TOPPADDING", (0, 0), (-1, -1), 2 * mm), ("BOTTOMPADDING", (0, 0), (-1, -1), 2 * mm),
            ("GRID", (0, 0), (-1, -1), 0.4, colors.lightgrey)]))

        story.append(table)
        if start + cells_per_page < len(image_files):
            story.append(PageBreak())
            story.append(Paragraph(title, title_style))

    doc.build(story)
    return len(image_files)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--best_poses_folder", required=True)
    parser.add_argument("--output_folder", default="DeepPK_Filtered_Results_Caco6_OneDILI",)
    args = parser.parse_args()

    input_csv = Path(args.input_csv).resolve()
    best_poses_folder = Path(args.best_poses_folder).resolve()
    output_folder = Path(args.output_folder).resolve()

    if not input_csv.is_file():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    if not best_poses_folder.is_dir():
        raise NotADirectoryError(f"Best poses folder not found: {best_poses_folder}")

    output_folder.mkdir(parents=True, exist_ok=True)
    best_pose_output = output_folder / "Extracted_Best_Poses"
    image_output = output_folder / "2D_Images"

    best_pose_output.mkdir(parents=True, exist_ok=True)
    image_output.mkdir(parents=True, exist_ok=True)

    df = normalize_dataframe_columns(pd.read_csv(input_csv))
    assessment = df.copy()
    rule_columns = []
    missing_columns = []

    for rule_name, (column_name, operator, threshold) in NUMERIC_RULES.items():
        column = find_column(df, column_name)
        if column is None:
            missing_columns.append(column_name)
            continue

        values = numeric_values(df[column])
        if operator == "<=":
            result = values <= threshold
        else:
            result = values >= threshold

        assessment[rule_name] = result.fillna(False)
        rule_columns.append(rule_name)

    for rule_name, (column_name, expected) in CATEGORICAL_RULES.items():
        column = find_column(df, column_name)
        if column is None:
            missing_columns.append(column_name)
            continue

        values = categorical_values(df[column])
        assessment[rule_name] = values.eq(expected)
        rule_columns.append(rule_name)

    if missing_columns:
        raise KeyError("Required columns are missing: " + ", ".join(sorted(set(missing_columns))))

    assessment["Criteria_Passed"] = assessment[rule_columns].sum(axis=1)
    assessment["Criteria_Total"] = len(rule_columns)
    assessment["ALL_CRITERIA_PASS"] = assessment[rule_columns].all(axis=1)
    assessment["Failed_Criteria"] = assessment.apply(
        lambda row: "; ".join(
            rule for rule in rule_columns if not bool(row[rule])
        ),
        axis=1,
    )

    strict_pass = assessment[assessment["ALL_CRITERIA_PASS"]].copy()
    failed = assessment[~assessment["ALL_CRITERIA_PASS"]].copy()

    assessment_path = output_folder / "DeepPK_filter_assessment_all_molecules.csv"
    strict_pass_path = output_folder / "DeepPK_strict_PASS_only.csv"
    failed_path = output_folder / "DeepPK_failed_compounds.csv"
    summary_path = output_folder / "DeepPK_filter_summary.txt"
    extraction_path = output_folder / "strict_pass_extraction_manifest.csv"
    missing_pose_path = output_folder / "missing_best_poses.txt"

    assessment.to_csv(assessment_path, index=False)
    strict_pass.to_csv(strict_pass_path, index=False)
    failed.to_csv(failed_path, index=False)

    extracted = []
    missing_poses = []
    image_failures = []

    for index, (_, row) in enumerate(strict_pass.iterrows(), start=1):
        number = str(row["Molecule_Number"]) if "Molecule_Number" in row.index else str(index)
        label = molecule_label(row, index)
        safe_label = safe_filename(label)
        identifiers = row_identifiers(row)

        pose_path = find_pose_file(str(best_poses_folder), identifiers)
        copied_pose = ""
        image_path = image_output / f"Molecule_{number}_{safe_label}.png"
        image_generated = False
        image_source = ""

        if pose_path is None:
            missing_poses.append(f"{number}\t{label}")
        else:
            destination = best_pose_output / os.path.basename(pose_path)
            counter = 1
            original_destination = destination

            while destination.exists():
                destination = best_pose_output / (
                    f"{original_destination.stem}_{counter}"
                    f"{original_destination.suffix}"
                )
                counter += 1

            shutil.copy2(pose_path, destination)
            copied_pose = str(destination)

            mol = mol_from_structure(pose_path)
            if mol is not None:
                try:
                    generate_2d_image(mol, image_path)
                    image_generated = True
                    image_source = "best_pose"
                except Exception:
                    image_generated = False

        if not image_generated:
            smiles = smiles_from_row(row)
            mol = Chem.MolFromSmiles(smiles) if smiles else None

            if mol is not None:
                try:
                    generate_2d_image(mol, image_path)
                    image_generated = True
                    image_source = "canonical_smiles_fallback"
                except Exception as exc:
                    image_failures.append(f"{number}\t{label}\t{exc}")

            else:
                image_failures.append(f"{number}\t{label}\tNo usable structure or SMILES")

        extracted.append(
            {
                "Molecule_Number": number,
                "Molecule_Label": label,
                "Best_Pose_File": copied_pose,
                "2D_Image": str(image_path) if image_generated else "",
                "2D_Image_Source": image_source,
                "Best_Pose_Found": bool(pose_path),
                "2D_Image_Generated": image_generated,
            }
        )

    manifest = pd.DataFrame(extracted)
    manifest.to_csv(extraction_path, index=False)

    with open(missing_pose_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(missing_poses))

    with open(output_folder / "image_generation_failures.txt", "w", encoding="utf-8") as handle:
        handle.write("\n".join(image_failures))

    pdf_path = output_folder / "Strict_Pass_Molecules_2D_Grid.pdf"
    pdf_image_count = create_2d_pdf(image_output, pdf_path, title="Strict-Pass Molecules - 2D Structures",)

    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write("Deep-PK filtering and strict-pass extraction summary\n")
        handle.write("====================================================\n\n")
        handle.write(f"Input CSV: {input_csv}\n")
        handle.write(f"Best poses folder: {best_poses_folder}\n")
        handle.write(f"Total molecules: {len(df)}\n")
        handle.write(f"Strict PASS: {len(strict_pass)}\n")
        handle.write(f"FAIL: {len(failed)}\n")
        handle.write(f"Criteria: {len(rule_columns)}\n")
        handle.write(f"Best poses extracted: {sum(bool(x['Best_Pose_Found']) for x in extracted)}\n")
        handle.write(f"Missing best poses: {len(missing_poses)}\n")
        handle.write(f"2D images generated: {sum(bool(x['2D_Image_Generated']) for x in extracted)}\n")
        handle.write(f"2D image failures: {len(image_failures)}\n")
        handle.write(f"2D image PDF: {pdf_path}\n")
        handle.write(f"Images included in PDF: {pdf_image_count}\n\n")

        handle.write("Active numeric criteria:\n")
        for _, (column, operator, threshold) in NUMERIC_RULES.items():
            handle.write(f"{column} {operator} {threshold}\n")

        handle.write("\nActive categorical criteria:\n")
        for _, (column, expected) in CATEGORICAL_RULES.items():
            handle.write(f"{column} = {expected}\n")

        handle.write("\nRemoved criteria:\n")
        handle.write("- Human Oral Bioavailability 20%\n")
        handle.write("- CYP 1A2, 2C19, 2C9, 2D6 and 3A4 inhibitor filters\n")
        handle.write("- BCRP\n")
        handle.write("- OATP1B1\n")
        handle.write("- OATP1B3\n")
        handle.write("- OCT2\n")
        handle.write("- Liver Injury II\n")
        handle.write("- DeepPK_ranked_by_criteria_passed.csv\n")

    print(f"Total molecules: {len(df)}")
    print(f"Strict PASS: {len(strict_pass)}")
    print(f"FAIL: {len(failed)}")
    print(f"Best poses extracted: {sum(bool(x['Best_Pose_Found']) for x in extracted)}")
    print(f"Missing best poses: {len(missing_poses)}")
    print(f"2D images generated: {sum(bool(x['2D_Image_Generated']) for x in extracted)}")
    print(f"2D PDF generated: {pdf_path}")
    print("DeepPK_ranked_by_criteria_passed.csv was not generated.")
    print(f"Output folder: {output_folder}")

if __name__ == "__main__":
    main()