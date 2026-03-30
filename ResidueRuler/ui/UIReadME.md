

# ResiRuler

This is ResiRuler's web interface for analyzing and visualizing residue-level structural changes in biomolecular models. It supports measuring distances between annotated residues, comparing distances across conformations, and visualizing movement vectors using ChimeraX.

---

## Features

- Extracts and compares residue distances from `.cif` structures
- Computes movement vectors between aligned conformations
- Generates ChimeraX scripts for structural visualization
- Plots difference and shift metrics for publication-quality graphics

---

### `run`: Distance Extraction Mode

![Run Page](../images/Run_Page.png)

Run mode will calculate pairwise distances between residues specified by an Excel Sheet and output of table of these calculations as well as a visualization script for chimeraX which will draw links between the residues and color those links according to  thresholds

**It requires**:
- Excel specifying the residues to calculate the distance between
- cif file of the atomic structure model
- JSON mapping specifying which label in the Excel sheet matches up with the chain IDs in the cif. The simplest way to figure out this mapping is through loading the cif into a molecular visualization software such as ChimeraX and visually inspecting the chain IDs to determine the proper mapping.
- JSON mapping specifying the thresholds to use for coloring the links in the visualization script

**Example Excel Format**:

![Excel Format Example](../images/ExcelEX.png)


After Inputting the proper files, it becomes easy to generate ChimeraX visualizations representing distances between important residues.

**Example Visualization**:

![Distance Links Example](../images/Link_EX.png)

---

### `compare`: Distance Difference Mode

![Compare Page](../images/Compare_Page.png)

Compare mode takes in two CSV outputs from Run mode and calculates how distances change between residues of interest between two structures. It also produces a bar plot representing the distance changes from one structure to the next 


**It Requires**:
- results CSV from Run mode on one structure
- results CSV from Run mode on another structure
- JSON mapping specifying which chain of the first structure corresponds to chains in the other structure if there are labeling differences

**Example Visualization**:

![Distance Difference Plot](../images/distance_diff_plot.png)

---

### `movement`: Movement Analysis Mode

![Movement Page](../images/Movement_Page.png)

Takes in two aligned structures/conformations and outputs the displacement for every residue with a mapping match from one conformation to the next 

**It Requires**:
- One aligned `.cif` file
- A second aligned `.cif` file
- JSON mapping for chains across models

**Example Visualizations**:

![Movement Shift Plot](../images/movement_plot.png)

![Movement Vectors](../images/movement_vectors.png)

![Chimera Movement Shift Plot](../images/ChimeraPlot_EX.png)

![Chimera Movement Vectors Plot](../images/Vector_EX.png)
---


## License

MIT License

---

## Author

Timothy Baker
