# Automated Pipeline
A typical SBVS workflow, run manually, spans six to eight disconnected tools: structure-repair software, a separate ligand-preparation toolkit, a docking engine, manual visual inspection of poses, an external website for ADMET prediction, and finally a spreadsheet for filtering.
Each transition between these tools is a place where file formats break, results get mismatched, and reproducibility becomes irretractable because no one can fully reconstruct, months later, exactly what was clicked and in what order.
This pipeline removes that fragmentation by chaining every stage into a single, traceable, file-to-file data flow.
A user supplies only two things: a target PDB ID (or structure file) and a list of PubChem CIDs and the pipeline carries that input through structure preparation, ligand preparation, automatic pocket detection, docking, energy filtering, interaction-aware pose selection, ADMET/PK prediction, and final rule-based filtering, without requiring the user to manually operate any of the underlying scientific software.

**What makes it genuinely “one-stop” in practice:
**

• Automatic structure sourcing and repair - no manual structure cleanup in a separate molecular viewer.

• Automatic pocket/grid definition from the co-crystallized ligand - no manual box-picking, which is one of the most error-prone manual steps in conventional docking workflows.

• Resumable, fault-tolerant docking - large libraries can be interrupted and restarted without losing completed work.

• A two-stage ranking system (energy, then interaction quality) that removes the need to manually eyeball every docked pose.

• A direct, automated bridge into Deep-PK ADMET prediction - no manual re-uploading of shortlisted molecules to an external website.

• A final, rule-based filter that outputs a ready-to-review shortlist, complete with 2D structures and a single consolidated PDF.
