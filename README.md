# Automated Pipeline
A typical SBVS workflow, run manually, spans six to eight disconnected tools: structure-repair software, a separate ligand-preparation toolkit, a docking engine, manual visual inspection of poses, an external website for ADMET prediction, and finally a spreadsheet for filtering.
Each transition between these tools is a place where file formats break, results get mismatched, and reproducibility becomes irretractable because no one can fully reconstruct, months later, exactly what was clicked and in what order.
This pipeline removes that fragmentation by chaining every stage into a single, traceable, file-to-file data flow.
A user supplies only two things: a target PDB ID and a list of PubChem CIDs and the pipeline carries that input through structure preparation, ligand preparation, automatic pocket detection, docking, energy filtering, interaction-aware pose selection, ADMET/PK prediction, and final rule-based filtering, without requiring the user to manually operate any of the underlying scientific software.

## What makes it genuinely “one-stop” in practice:

• Automatic structure sourcing and repair - no manual structure cleanup in a separate molecular viewer.

• Automatic pocket/grid definition from the co-crystallized ligand - no manual box-picking, which is one of the most error-prone manual steps in conventional docking workflows.

• Resumable, fault-tolerant docking - large libraries can be interrupted and restarted without losing completed work.

• A two-stage ranking system (energy, then interaction quality) that removes the need to manually eyeball every docked pose.

• A direct, automated bridge into Deep-PK ADMET prediction - no manual re-uploading of shortlisted molecules to an external website.

• A final, rule-based filter that outputs a ready-to-review shortlist, complete with 2D structures and a single consolidated PDF.

# Pipeline Workflow at a glance

<img width="818" height="983" alt="image" src="https://github.com/user-attachments/assets/f4c0351f-4286-442c-b6bd-c841f3c6ea4d" />

# System Requirements, Environment Setup and Execution Protocol

Windows/PowerShell environment and the ordered commands used to execute the workflow.

# Operating-system and hardware requirements

•	Windows with PowerShell and a working Conda installation (Anaconda or Miniconda).

•	Internet access for PDB/PDB-REDO, PubChem and Deep-PK website access.

•	AutoDock Vina installed and accessible from PowerShell.

•	MGLTools/AutoDockTools available for receptor PDBQT preparation.

•	A supported Chrome browser for browser-driven Deep-PK submissions.

•	Adequate disk space for ligand libraries, docking outputs, logs, downloaded result tables, images and PDF reports.

# Software and Package Requirements
1. Python 3.11 + Conda - Python runtime and environment management
2. OpenMM and PDBFixer - Protein structure processing and repair
3. RDKit - Ligand handling, SMILES conversion and 2D drawing
4. Meeko - Ligand PDBQT preparation
5. Gemmi - Structure-file parsing where used
6. AutoDock Vina - Molecular docking executable
7. MGLTools/AutoDockTools - Receptor PDBQT preparation
8. Playwright - Browser automation; install the library imported by the local DeepPK_results.py
9. pandas, Pillow, ReportLab - Tabular processing, image handling and PDF generation where used

# Create Main Conda Environment
The following commands creates the environment and install all the packages required for the entire workflow.

```
conda create -n docking python=3.11
conda activate docking
conda install -c conda-forge openmm pdbfixer rdkit meeko gemmi pandas pillow reportlab
python -m pip install playwright
python -m playwright install chromium
```
# Optional separate RDKit drawing environment
Use this separate environment only if you want to isolate/test RDKit drawing dependencies.

```
conda create -n rdkit_draw -c conda-forge python=3.11 rdkit pandas pillow
conda activate rdkit_draw
python -c "from rdkit import Chem; from rdkit.Chem import Draw; print('RDKit 2D drawing OK')"
```
# Steps for Running the Entire Protocol

## Step 1: Activate the environment and enter the project folder
Check the prompt shows the docking environment and the correct working directory. 

```
conda activate docking
cd "E:\Automation\PDB_Redo"
```
## Step 2: Prepare the target protein
Enter the correct PDB ID. A more refined structure can also be downloaded from PDBRedo if the given PDB ID is registered in the PDBRedo database. Make sure that the correct paths to MGL_PYTHON and PREPARE_RECEPTOR in the code are mentioned. Make sure that correct path to the directory where MGLTools have been installed. 

```
python .\pdb_redo.py
```
After running the script, verify the cleaned/selected PDB and receptor PDBQT outputs exist.

## Step 3: Download and prepare the ligand library
A summary file needs to be downloaded from PubChem database containing CIDS.

```
python .\PubChem.py
```
Verify SDF files in Molecules foler and prepared ligand PDB/PDBQT folders are created with the PDB and PDBQT files resectively.

## Step 4: Generate the grid box and target-residue list
If a bound ligand is present in the structure is present, a gridbox will be defined and interacting residues (residues.txt) file will be generated.

```
python .\Grid_gen.py 2VHW_PDBREDO_selected.pdb
```
Verify gridbox.txt, vina_config.txt, interacting_residues.csv and residues.txt.

## Step 5: Run AutoDock Vina
Place the vina.exe file in the working directory. 

```
python .\Vina_run.py
```
Verify Docking_Results contains per-ligand out.pdbqt and log.txt, and inspect docked/undocked summary CSVs.

## Step 6: Extract binding energies and apply the cutoff
Enter the cuttoff binding energy without units. Binding energy of the given value and lower will be reported

```
python .\BE_data.py
```
Verify binding_energies.csv, cutoff_binding_energies.csv and Cutoff_Results.

## Step 7: Select best poses and create the top-screened list
Based on the interacting residues given in the residues.txt file are checked for the all the 9 poses of the docked molecules. If the molecules show consistent interactions to those mentioned in the residues.txt file are mentioned and reported with the best pose and type on interactions in the interactions_summary.csv file.

```
python .\Best_pose.py
```
It generates Best_Poses, interaction summaries and top10_screened.csv.

## Step 8: Run Deep-PK on the selected molecules
DeepPK_results.py uses top10_screened.csv to determine which molecules should enter the ADMET stage. It then reads the original SDF folder and matches selected records to their original structures. This is important because the original SDF retains molecular connectivity and stereochemical information more reliably than reconstructing a molecule from a docking PDB. Each selected SDF molecule is converted to canonical isomeric SMILES using RDKit. A SMILES input file is generated for Deep-PK. Metadata records molecule number, SDF file, canonical SMILES and source information. selected_molecules.json preserves the mapping used during the run.

```
python .\DeepPK_results.py --sdf_folder "E:\Automation\PDB_Redo\molecules" --top10_csv "E:\Automation\PDB_Redo\top10_screened.csv" --output_folder "E:\Automation\PDB_Redo\DeepPK"
```
Verify ADMET submission completes and the main results CSV is saved in DeepPK.

## Step 9: Filter Deep-PK results
Based on the drug-likeness and ADMET properties the molecules were screened. The parameters are selective and can be changed upon user requirements directly in the code.

```
conda activate rdkit_draw
python -c "from rdkit import Chem; from rdkit.Chem import Draw; print('RDKit 2D drawing OK')"
python .\DeepPK_filter_molecules.py --input_csv "E:\Automation\PDB_Redo\DeepPK\DeepPK_all_properties.csv" --output_folder "E:\Automation\PDB_Redo\Filtered_mols"
```
If --input_csv expects a file, pass DeepPK_all_properties.csv rather than the DeepPK folder.









