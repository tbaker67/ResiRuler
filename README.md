<div align="center">

# ResiRuler

**A web-based tool for residue-level structural comparison of biomolecular models**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](#license)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](#requirements)
[![Streamlit](https://img.shields.io/badge/UI-Streamlit-FF4B4B.svg)](#running-the-ui)

</div>

---

ResiRuler analyzes residue-level structural changes across biomolecular models, aligning structures and quantifying how residues move relative to one another. Results can be exported and viewed directly in **ChimeraX** and **PyMOL**.

It Features:
- **Alignment** of structures using [USalign's](https://www.nature.com/articles/s41592-022-01585-1) MMAlign.
- **Automatic chain mapping** between reference and comparison structures, tunable by sequence identity, RMSD penalty, and minimum percent identity.
- **Pairwise distance comparison**, showing how distances and contacts between residue pairs change across structures.
- **Per-residue displacement**, visualizing the direction and magnitude of residue movement between structures.
- **Interactive, downloadable plots** — contact maps and distance-difference maps all exportable as SVGs.
- **Ready-made visualization scripts** for ChimeraX and PyMOL.

---

## Requirements

| Requirement | Details |
|---|---|
| OS | Linux or macOS |
| Python | 3.11 |
| Package manager | [Conda](https://www.anaconda.com/docs/getting-started/getting-started) |

---

## Installation

```bash
git clone https://github.com/tbaker67/ResiRuler.git
cd ResiRuler
conda env create -f environment.yml
conda activate resiruler
```
---

## Running the UI

ResiRuler's web interface is managed through **Streamlit**.

1. Activate the conda environment:

   ```bash
   conda activate resiruler
   ```

2. Launch the UI:

   ```bash
   streamlit run ResidueRuler/ui/app.py
   ```

A successful launch prints a message like the one below. The **Network URL** lets others on the same local network reach their own independent instance of the interface.

![Successful Launch](ResidueRuler/images/Streamlit_success.png)

---

## The "Align" Tab

Align two related structures with MMAlign before comparing them, or skip this tab entirely if you already have prealigned structures.

<details>
<summary>Show steps</summary>

1. **Upload your files.**

   ![Align Page 1](ResidueRuler/images/align/Align_Page1.png)

2. **Select the chains** to include in the alignment. Selecting fewer chains produces a faster alignment.

   ![Align Page 2](ResidueRuler/images/align/Align_Page2.png)

3. **Preview and download** the aligned result.

   ![Align Page 3](ResidueRuler/images/align/Align_Page3.png)

4. **(Optional) Clean the alignment** to drop unmatched chains or residues, and adjust the parameters used to determine matched chains and residues.

   ![Align Page 4](ResidueRuler/images/align/Align_Page4.png)

</details>

---

## The "Pairwise Distance Difference" Tab

See how distances between pairs of residues change across structures, represented as both raw distance changes and changes in "contacts" based on thresholds you define.

<details>
<summary>Show steps</summary>

1. **Upload your aligned structures**: one reference structure, plus as many comparison structures as you want.

   ![Pairwise Distance Difference Tab 1](ResidueRuler/images/compare/Compare_Page1.png)

2. **Define chain correspondence** between the reference and comparison structures. Any chain you don't explicitly match is handled by the automatic chain-mapping algorithm, which you can tune via settings for sequence alignment, RMSD penalty, and minimum percent identity required for a valid match.

   > Hover over any question-mark icon in the UI for more detail on these settings.

   ![Pairwise Distance Difference Tab 2](ResidueRuler/images/compare/Compare_Page2.png)

3. **Set distance thresholds** for what gets displayed in the output.

   ![Pairwise Distance Difference Tab 3](ResidueRuler/images/compare/Compare_Page3.png)

4. **Select chains to display.** A panel with contact maps, distance-difference maps, and a data table appears. All plots are interactive, resizable, and downloadable as SVGs.

   ![Pairwise Distance Difference Tab 4](ResidueRuler/images/compare/Compare_Page4.png)

</details>

---

## The "Per-Residue Displacement" Tab

Calculate the displacement between residues in the reference structure and their corresponding counterparts in a comparison structure, and visualize the direction and magnitude of that movement.

<details>
<summary>Show steps</summary>

1. **Upload structures and define chain correspondence**, just as in the Pairwise Distance Difference tab.

   ![Per-Residue Displacement Tab 1](ResidueRuler/images/compare/Compare_Page1.png)
   ![Per-Residue Displacement Tab 2](ResidueRuler/images/compare/Compare_Page2.png)

2. **Run the analysis.** Click "Analyze Displacement," then choose a custom color scale and adjust the vector fidelity slider to control how many displacement vectors are shown.

   ![Per-Residue Displacement Tab 3](ResidueRuler/images/movement/Displacement_Page1.png)

3. **Review the results.** In-app visualizations and a full data table are displayed.

   ![Per-Residue Displacement Tab 4](ResidueRuler/images/movement/Displacement_Page2.png)

4. **Download your results** as a zipped folder containing the data, visualization scripts, and models.

   ![Per-Residue Displacement Tab 5](ResidueRuler/images/movement/Displacement_Page3.png)
   ![Per-Residue Displacement Tab 6](ResidueRuler/images/movement/Displacement_Page4.png)

</details>

### Downloaded Script Reference

Each displacement download includes ready-made scripts:

| File | Purpose |
|---|---|
| `coloring_script.pml` | Colors all models in **PyMOL** |
| `chimera_coloring_script.cxc` | Colors all models in **ChimeraX** |
| `pml_arrows/run_all_arrows.pml` | Displays all displacement arrows in **PyMOL** |
| `bild/open_all_bilds.cxc` | Displays all displacement arrows in **ChimeraX** |

> Models are named `{model being shown}-{model the colors show displacements against}` within PyMol and ChimeraX.

**Example visualization:**

![Example Per-Residue Displacement Tab Visualization](ResidueRuler/images/movement/Chimera_Movemement_Map.png)

### Recoloring in ChimeraX

Recolor the displacement gradient inside ChimeraX via **Tools → Depiction → Render/Select by Attribute**.

![Per-Residue Displacement Tab Edit 1](ResidueRuler/images/movement/Movement_Edit1.png)

Select the `distance` attribute in the window that appears to threshold and recolor with a different palette.

![Per-Residue Displacement Tab Edit 2](ResidueRuler/images/movement/Movement_Edit2.png)

---

## Remote Access to ResiRuler on an HPC Cluster

> Example port: `8501` · Example user: `user` — substitute your own port and username below.

<details>
<summary>Show steps</summary>

1. **Start an interactive session:**

   ```bash
   srun --partition=sb-gpu --nodes=1 --ntasks=1 --cpus-per-task=12 \
        --mem-per-cpu=3g --gpus=2080ti:2 --time 4:00:00 --pty /bin/bash
   ```

2. **Check which GPU node you've been assigned** (e.g. `gpucomp-01`), then activate the environment and launch the app:

   ```bash
   hostname
   conda activate resiruler
   cd path/to/ResiRuler/ui/
   streamlit run app.py
   ```

3. **In a new terminal window, forward the port:**

   ```bash
   ssh -L 8501:localhost:8501 user@login.edu
   ssh -N -L 8501:localhost:8501 gpucomp-01
   ```

   Enter your password when prompted — the terminal will appear to hang while the tunnel is active, which is expected.

4. **In a browser on your local machine**, navigate to:

   ```
   localhost:8501
   ```

</details>

---

## License

Distributed under the [MIT License](LICENSE).

---

## Author

Timothy Baker, Wilhelm Salmen