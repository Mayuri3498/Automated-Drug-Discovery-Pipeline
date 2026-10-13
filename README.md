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


