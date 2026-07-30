"""Utility functions and UI widgets for the ResiRuler Streamlit interface."""

import gc
import io
import os
import random
import tempfile
import zipfile
from contextlib import contextmanager
from io import BytesIO

import numpy as np
import pandas as pd
import streamlit as st
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from src.resiruler.core.auto_alignment import EnsembleMapper, StructureMapper
from src.resiruler.core.structure_parsing import (
    extract_res_from_chain,
    load_structure,
)

DIVERGING_COLORSCALES = [
    "RdBu_r",
    "RdBu",
    "PiYG",
    "PRGn",
    "Spectral",
    "Picnic",
    "Portland",
    "Balance",
]

SEQUENTIAL_COLORSCALES = [
    "Viridis",
    "Plasma",
    "Cividis",
    "Greens",
    "Reds",
    "Blues",
    "Oranges",
    "YlOrRd",
]

def reset_downstream(*keys):
    """
    Drop stale large objects (matrices, structures, figures, CIF strings, etc.)
    from session_state and force a garbage-collection pass.
    """
    freed_any = False
    for key in keys:
        if key in st.session_state:
            del st.session_state[key]
            freed_any = True
    if freed_any:
        gc.collect()


@contextmanager
def struct_to_temp_cif(structure):
    """
    Create a temporary CIF file for a Bio.PDB structure.
    Automatically deletes the file after use.
    """
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".cif")
    try:
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(tmp.name)
        tmp.close()
        yield tmp.name
    finally:
        if os.path.exists(tmp.name):
            os.unlink(tmp.name)


def save_temp_file(uploaded_file):
    with tempfile.NamedTemporaryFile(
        delete=False, suffix=os.path.splitext(uploaded_file.name)[-1]
    ) as tmp_file:
        tmp_file.write(uploaded_file.getbuffer())
        return tmp_file.name


def create_downloadable_zip(files_dict):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as z:
        for filename, content in files_dict.items():
            z.writestr(filename, content)
    buffer.seek(0)
    return buffer


def create_ensemble_mapper(
    ref_structure,
    tgt_structures,
    chain_mappings,
    threshold,
    protein_aligner,
    nucleotide_aligner,
    rmsd_scaling_factor
):
    ensemble_mapper = EnsembleMapper(
        ref_structure,
        protein_aligner=protein_aligner,
        nucleotide_aligner=nucleotide_aligner,
    )

    for structure_name, tgt_structure in tgt_structures.items():
        ensemble_mapper.add_structure(
            structure_name,
            tgt_structure,
            threshold,
            chain_mappings[structure_name],
            rmsd_scaling_factor
        )

    return ensemble_mapper


def chain_selector_ui(
    structure,
    label="Select Chains to Visualize (Chains with All will be deselected)",
    default_all=True,
    key_prefix="chain_selector",
):
    """
    Display a Streamlit multiselect box for choosing chains from a structure or structure mapper.
    If "All" is selected along with individual chains, returns all chains EXCEPT those individual chains.
    As such, structure should only be a BioPython structure object, or a StructureMapper object from autoalignments.py
    """
    key = key_prefix + "chain_selector"
    if structure is None:
        st.warning("No structure loaded. Please upload a file first.")
        return None

    chain_options = None
    if isinstance(structure, StructureMapper):
        chain_options = sorted(structure.matched_ref_chains)
    else:
        chain_options = sorted({chain.id for model in structure for chain in model})
    chain_options = sorted(chain_options)

    selected = st.multiselect(
        label,
        chain_options,
        default=[],
        key=f"{key_prefix}_chain_order_selector",
        help=(
            "Order here controls the order chains appear along the heatmap axes. "
            "Remove and re-add a chain to move it to the end."
        ),
    )
    return selected


def load_structure_if_new(cif_file, name_key, struct_key):
    if cif_file is None:
        return None

    if name_key not in st.session_state or st.session_state[name_key] != cif_file.name:
        reset_downstream(
            struct_key,
            "mapper",
            "ref_dm",
            "tgt_dms_dict",
            "compare_dms_dict",
            "small_ref_dm",
            "small_tgt_dms_dict",
            "small_compare_dms_dict",
            "ref_fig",
            "tgt_figs_dict",
            "compare_figs_dict",
            "displacement_dfs",
        )

        cif_path = save_temp_file(cif_file)
        structure = load_structure(cif_path)

        st.session_state[struct_key] = structure
        st.session_state[name_key] = cif_file.name
    else:
        # same file
        structure = st.session_state.get(struct_key)

    return structure


def load_structures_if_new(cif_files, name_key_prefix, struct_key_prefix):
    """
    Load one or more structures, caching them in session_state.
    Returns a dict mapping file name -> structure object.
    """
    if not cif_files:
        return {}

    structures = {}
    any_changed = False
    for i, cif_file in enumerate(cif_files):
        name_key = f"{name_key_prefix}_{i}"
        struct_key = f"{struct_key_prefix}_{i}"

        if (
            name_key not in st.session_state
            or st.session_state[name_key] != cif_file.name
        ):
            any_changed = True
            # free the structure this index used to hold before overwriting it
            st.session_state.pop(struct_key, None)

            cif_path = save_temp_file(cif_file)
            structure = load_structure(cif_path)

            st.session_state[struct_key] = structure
            st.session_state[name_key] = cif_file.name
        else:
            structure = st.session_state.get(struct_key)

        # remove .cif for names
        structures[cif_file.name[:-4]] = structure

    # if fewer files are uploaded than before, drop the now-unused cached
    i = len(cif_files)
    while f"{name_key_prefix}_{i}" in st.session_state:
        st.session_state.pop(f"{name_key_prefix}_{i}", None)
        st.session_state.pop(f"{struct_key_prefix}_{i}", None)
        any_changed = True
        i += 1

    if any_changed:
        # target set changed so reset all things downstream
        reset_downstream(
            "mapper",
            "ref_dm",
            "tgt_dms_dict",
            "compare_dms_dict",
            "small_ref_dm",
            "small_tgt_dms_dict",
            "small_compare_dms_dict",
            "ref_fig",
            "tgt_figs_dict",
            "compare_figs_dict",
            "displacement_dfs",
        )

    return structures


def chain_mapping_input(ref_chains, tgt_chains, default=None, key="chain_mapping"):
    """
    Editable table for mapping reference → target chains.
    """
    st.caption(
        "Map reference chains → target chains (leave all blank for auto mapping)"
    )

    # Build default DataFrame
    with st.expander("Explicit Chain Mapping Options"):
        df = pd.DataFrame(
            [(ref, "", "") for ref in ref_chains],
            columns=["Reference Chain", "Target Chain", "Chain Type"],
        )

        # Editable table
        edited_df = st.data_editor(
            df,
            num_rows="fixed",
            key=key,
            column_config={
                "Reference Chain": st.column_config.TextColumn(
                    "Reference Chain", disabled=True
                ),
                "Target Chain": st.column_config.SelectboxColumn(
                    "Target Chain",
                    options=[None] + list(tgt_chains),  # allow blank
                    required=False,
                ),
                "Chain Type": st.column_config.SelectboxColumn(
                    "Chain Type", options=[None] + ["protein", "dna", "rna"]
                ),
            },
            use_container_width=True,
            hide_index=True,
        )

        # Convert to dict
        mapping = {
            ref: (
                None if pd.isna(tgt) or tgt == "" else tgt,
                (None if pd.isna(chain_type) or chain_type == "" else chain_type),
            )
            for ref, tgt, chain_type in zip(
                edited_df["Reference Chain"],
                edited_df["Target Chain"],
                edited_df["Chain Type"],
            )
        }
        return mapping


def get_chain_mappings_for_targets(
    structures_dict, ref_chains, default_mapping=None, key="mapping"
):
    """
    Build chain mapping table for each target structure.
    """
    mappings = {}
    for filename, tgt_structure in structures_dict.items():
        st.subheader(f"Chain mapping for target: {filename}")

        tgt_chains = [chain.id for chain in tgt_structure.get_chains()]
        mapping = chain_mapping_input(
            ref_chains,
            tgt_chains,
            default=default_mapping,
            key=f"{key}_{filename}",
        )
        mappings[filename] = mapping

    return mappings


def get_threshold(label, default, key):
    threshold_input = st.text_input(label, value=default, key=key,
                                    help="Minimum Percent Identity required for any two chains to be matched with one another"
                                    )
    try:
        return float(threshold_input)
    except ValueError:
        st.error("Please enter a valid number.")
        return None
    


def create_downloadable_zip_grouped(file_groups):
    """
    Create a ZIP buffer from multiple file groups.
    file_groups = dict where
        key   = folder name or None for root
        value = dict mapping {filename: content}

    Special case:
        key=None → files go in the root of the ZIP
        key="csv"/"bild"/"models" → files go in that subfolder
    Supports str, bytes, StringIO, BytesIO as content.
    """
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as z:
        for folder, files in file_groups.items():
            for filename, content in files.items():
                # Convert buffers to bytes
                if hasattr(content, "getvalue"):
                    content = content.getvalue()
                if isinstance(content, str):
                    content = content.encode("utf-8")
                # Determine path inside ZIP
                arcname = f"{folder}/{filename}" if folder else filename
                z.writestr(arcname, content)
    buffer.seek(0)
    return buffer


PROTEIN_SUBSTITUTION_MATRICES = [
    "BLOSUM62",
    "BLOSUM45",
    "BLOSUM50",
    "BLOSUM80",
    "BLOSUM90",
    "BLASTP",
    "DAYHOFF",
    "FENG",
    "GENETIC",
    "GONNET1992",
    "JOHNSON",
    "JONES",
    "LEVIN",
    "MCLACHLAN",
    "MDM78",
    "BENNER22",
    "BENNER6",
    "BENNER74",
    "PAM250",
    "PAM30",
    "PAM70",
    "RAO",
    "RISLER",
]

NUCLEOTIDE_SUBSTITUTION_MATRICES = ["MEGABLAST", "BLASTN", "NUC.4.4", "HOXD70"]


def aligner_ui(protein, key_prefix="aligner"):
    """
    Function to allow user to specify a biopython Pairwise Aligner Object
    """
    with st.container():
        # Alignment mode and matrix
        col1, col2 = st.columns(2)
        with col1:
            mode = st.selectbox(
                "Alignment Mode", ["global", "local"], key=f"{key_prefix}_mode",
                help="Whether to use the whole alignment (global), "
                "or the best matching sub-region (local)."
            )
        with col2:
            if protein:
                matrix_name = st.selectbox(
                    "Substitution Matrix",
                    PROTEIN_SUBSTITUTION_MATRICES,
                    key=f"{key_prefix}_matrix",
                    help="The substitution matrix that will be used during alignment."
                )
            else:
                matrix_name = st.selectbox(
                    "Substitution Matrix",
                    NUCLEOTIDE_SUBSTITUTION_MATRICES,
                    key=f"{key_prefix}_matrix",
                    help="The substitution matrix that will be used during alignment."
                )
        # Gap scores in one row
        col3, col4 = st.columns(2)
        with col3:
            open_gap_score = st.number_input(
                "Gap Opening Penalty",
                value=-10,
                step=1,
                key=f"{key_prefix}_open_gap",
                help="Value added to the alignment score when a gap is opened." \
                " A more negative value means that gaps are less likely to be opened in the alignment."

            )
        with col4:
            extend_gap_score = st.number_input(
                "Gap Extension Penalty",
                value=-1,
                step=1,
                key=f"{key_prefix}_extend_gap",
                help="Value subtracted from the alignment score when an existing gap is made longer." \
                " A more negative value means that gaps are less likely to be extended in the alignment."
            )
        
        # optional end-gap penalties
        with st.expander("Advanced: End-gap penalties"):
            use_end_gaps = st.checkbox(
                "Enable End-Gap Penalties?", key=f"{key_prefix}_use_end_gaps",
                help="Whether or not to apply penalties to gaps that occur before (left) or after (right) a sequence."
            )
            left_open_gap_score = None
            right_open_gap_score = None
            if use_end_gaps:
                col5, col6 = st.columns(2)
                with col5:
                    left_open_gap_score = st.number_input(
                        "Left End Gap Penalty",
                        value=0,
                        step=1,
                        key=f"{key_prefix}_left_end_gap",
                        help="Value added to the alignment score for any gaps that occur at the beginning of the sequence." \
                        " A more negative values means that gaps at the beginning of the sequence are less likely."
                    )
                with col6:
                    right_open_gap_score = st.number_input(
                        "Right End Gap Penalty",
                        value=0,
                        step=1,
                        key=f"{key_prefix}_right_end_gap",
                        help="Value added to the alignment score for any gaps that occur at the end of the sequence." \
                        " A more negative values means that gaps at the end of the sequence are less likely."
                    )

    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = substitution_matrices.load(matrix_name)
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    if use_end_gaps:
        aligner.open_left_gap_score = left_open_gap_score
        aligner.open_right_gap_score = right_open_gap_score

    return aligner


def full_aligner_ui(key):
    with st.container():
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Protein Pairwise Aligner Settings")
            protein_aligner = aligner_ui(
                protein=True, key_prefix=f"{key} protein compare aligner"
            )
        with col2:
            st.subheader("Nucleotide Pairwise Aligner Settings")
            nucleotide_aligner = aligner_ui(
                protein=False, key_prefix=f"{key} nucleotide compare aligner"
            )

        rmsd_score_scaling_factor = st.number_input(
            "Scaling Factor For RMSD Penalty",
            value = float(1e-6),
            key=f"{key}_rmsd_scaling_factor",
            min_value=0.0,
            step=1e-6,
            format="%.8f",
            help= "This value is used to adjust how much to " \
                "penalize rmsd between chains during the chain matching process. " \
                "A higher value means that chains that are further apart will be more heavily penalized"
        )
    return protein_aligner, nucleotide_aligner, rmsd_score_scaling_factor


def show_alignments(ensemble_mapper, key="alignment_select"):
    """
    shows the alignment present in all structure mappers of a selected chiain in the reference
    """
    selected_chain = st.selectbox(
        "Select Chain to View Alignment of",
        options=[chain.id for chain in ensemble_mapper.ref_structure.get_chains()],
        key=key,
    )

    for structure_name, mapper in ensemble_mapper.structure_mappings.items():
        chain_mapper = mapper.chain_mappings.get(selected_chain, None)
        if chain_mapper is None or chain_mapper.alignment is None:
            continue

        st.subheader(
            f"Alignment for structure: {structure_name}, chain: {selected_chain}"
        )

        show_chain_alignment(chain_mapper)


def show_chain_alignment(chain_mapper):
    """
    Creates html for an alignment viewer, with highlighted colors indicating mismatches(red), gaps(grey), and matches(green)
    additionally, by hovering over a residue in the sequence the residue number will be shown
    """
    if chain_mapper is None or chain_mapper.alignment is None:
        st.warning("No Alignment Available")
        return

    seqs = [seq for seq in chain_mapper.alignment]
    if len(seqs) < 2:
        st.warning("Alignment has less than 2 sequences")
        return

    ref_seq = chain_mapper.aligned_ref_seq
    tgt_seq = chain_mapper.aligned_tgt_seq
    ref_res = extract_res_from_chain(chain_mapper.ref_chain)
    tgt_res = extract_res_from_chain(chain_mapper.tgt_chain)

    idx_ref, idx_tgt = 0, 0

    ref_line, tgt_line = "", ""

    for a, b in zip(ref_seq, tgt_seq):
        # color based on matches
        color = (
            "lightgray"
            if a == "-" or b == "-"
            else "lightgreen"
            if a == b
            else "lightcoral"
        )

        # reference residue
        if a != "-":
            res = ref_res[idx_ref]
            # includes insertion code if present
            ref_num = f"{res.id[1]}{res.id[2].strip()}"
            idx_ref += 1
            ref_line += (
                f"<span class='residue' style='background-color:{color}'>"
                f"{a}<span class='tooltip'>{ref_num}</span></span>"
            )
        else:
            ref_line += (
                f"<span class='residue' style='background-color:{color}'>{a}</span>"
            )

        # target residue
        if b != "-":
            res = tgt_res[idx_tgt]
            tgt_num = f"{res.id[1]}{res.id[2].strip()}"
            idx_tgt += 1
            tgt_line += (
                f"<span class='residue' style='background-color:{color}'>"
                f"{b}<span class='tooltip'>{tgt_num}</span></span>"
            )
        else:
            tgt_line += (
                f"<span class='residue' style='background-color:{color}'>{b}</span>"
            )

    # Creates an html object that can be viewed in streamlit
    html_content = f"""
    <html>
    <head>
    <style>
    .residue {{
        position: relative;
        display: inline-block;
        width: 1.3ch;  /* fixed width per residue */
        font-family: monospace;
        text-align: center;
        cursor: pointer;
    }}
    .tooltip {{
        display: none;
        position: absolute;
        top: -0.8em;
        left: 50%;
        transform: translateX(-50%);
        background: #eee;
        border: 1px solid #999;
        font-size: 10px;
        padding: 1px 3px;
        border-radius: 3px;
        white-space: nowrap;
        z-index: 10;
    }}
    .residue:hover .tooltip {{
        display: block;
    }}
    .sequence {{
        white-space: pre;
    }}
    .label {{
        display: inline-block;
        width: 80px;   /* fixed width for both labels */
        text-align: right;
        font-weight: bold;
    }}
    </style>
    </head>
    <body>
        <div class='sequence'><span class='label'>Ref ({chain_mapper.ref_chain.id}):</span> {ref_line}</div>
        <div class='sequence'><span class='label'>Tgt ({chain_mapper.tgt_chain.id}):</span> {tgt_line}</div>
    </body>
    </html>
    """

    st.components.v1.html(html_content, height=50, scrolling=True)


def get_measurement_mode(key="measurement_mode"):
    st.markdown("Choose Where to Measure Distances From")
    protein_option_to_mode = {
        "C-α": "CA",
        "C-β": "CB",
        "Side-Chain Centroid (Heavy Atoms Excluding CA)": "SC",
    }
    nucleic_option_to_mode = {"C1'": "C1'", "Centroid ()": "NC"}

    col1, col2 = st.columns(2)
    with col1:
        selected_protein_mode = st.selectbox(
            "Choose Which Atoms To Measure Proteins From",
            protein_option_to_mode.keys(),
            key=f"{key}_protein",
        )
    with col2:
        selected_nucleic_mode = st.selectbox(
            "Choose Which Atoms To Measure Nucleic Acids From",
            nucleic_option_to_mode.keys(),
            key=f"{key}_nucleic",
        )

    return (
        protein_option_to_mode[selected_protein_mode],
        nucleic_option_to_mode[selected_nucleic_mode],
    )


def filter_df_by_chains(df, selected_chains):
    """Filter dataframe to only include residues from selected chains in the reference structure."""
    if df.empty or not selected_chains:
        return df
    chain1 = df["ChainID_Resnum1"].str.split("-").str[0]
    mask = chain1.isin(selected_chains)
    return df[mask].reset_index(drop=True)


def distance_threshold_ui(key_prefix="distance_threshold",
                          min=float(0),
                          max=float(100)
                        ):
    """
    Function to allow user input for distance thresholding with a double-ended slider.
    """
    enable_threshold = st.checkbox(
        "Enable Contact Filtering",
        value=False,
        help="Filter displayed contacts by distance range. Shows shared contacts, gained, and lost.",
        key=f"{key_prefix}_enable",
    )

    if not enable_threshold:
        return False, None, None

    st.markdown("**Distance Range (Å)**")

    col_low, col_high = st.columns(2)
    with col_low:
        lower_threshold = st.number_input(
            "Lower (Å)",
            step=0.5,
            key=f"{key_prefix}_lower_num",
        )
    with col_high:
        upper_threshold = st.number_input(
            "Upper (Å)",
            value=float(15),
            step=float(0.5),
            key=f"{key_prefix}_upper_num",
        )

    if lower_threshold >= upper_threshold:
        st.warning("Lower threshold must be less than upper threshold.")

    st.caption(
        f"Showing Contacts With Distances Between {lower_threshold:.1f}Å and {upper_threshold:.1f}Å"
    )

    return enable_threshold, lower_threshold, upper_threshold


def color_scale_ui(purpose="diverging", key_prefix="color_scale", default_colorscale=None, min=float(0), max=float(100)):
    """
    Widget for letting the user customize a plot's color scale.
    """
    palette_options = DIVERGING_COLORSCALES if purpose == "diverging" else SEQUENTIAL_COLORSCALES
    default = default_colorscale or palette_options[0]
    label = "Difference (ΔDistance) Maps" if purpose == "diverging" else "Distance Maps"

    with st.expander(f"Color Scale Settings — {label}", expanded=False):
        colorscale = st.selectbox(
            "Colorscale",
            palette_options,
            index=palette_options.index(default) if default in palette_options else 0,
            key=f"{key_prefix}_colorscale",
        )

        vmin, vmax = None, None
        col1, col2 = st.columns(2)
        with col1:
            vmin = st.number_input(
                    "Min Value For Color Scale (Å)",
                    min_value=min,
                    value=min,
                    step=float(0.5),
                    key=f"{key_prefix}_vmin",
                )
        with col2:
            vmax = st.number_input(
                "Max Value For Color Scale (Å)",
                min_value=min,
                value=max,
                step=float(0.5),
                key=f"{key_prefix}_vmax",
            )

    return {"colorscale": colorscale, "vmin": vmin, "vmax": vmax}


def get_chain_pair_options(index_map):
    """Generate list of chain pair options from index_map."""
    chains = sorted({chain for chain, _ in index_map})
    pairs = []
    for i, c1 in enumerate(chains):
        for c2 in chains[i:]:  # Include self-pairs and upper triangle
            if c1 == c2:
                pairs.append((c1, c2, f"{c1} (intra-chain)"))
            else:
                pairs.append((c1, c2, f"{c1} × {c2}"))
    return chains, pairs


def add_svg_download(fig, filename):
    svg_buffer = BytesIO()
    fig.write_image(svg_buffer, format="svg")

    st.download_button(
        label="Download SVG",
        data=svg_buffer.getvalue(),
        file_name=filename,
        mime="image/svg+xml",
        key=f"download_{random.random()}",
    )


def get_available_chains(structure_path):
    """
    Gets available chains from a structure path.
    """
    parser = MMCIFParser(QUIET=True)
    structure_id = os.path.splitext(os.path.basename(structure_path))[0]
    structure = parser.get_structure(structure_id, structure_path)
    return sorted([chain.id for model in structure for chain in model])


def select_chains_for_alignment_ui(structure_path_1, structure_path_2):
    available_chains_1 = sorted(get_available_chains(structure_path_1))
    available_chains_2 = sorted(get_available_chains(structure_path_2))

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Chains to Consider in Reference**")
        selected_chains_1 = st.multiselect(
            "Chains to Align",
            available_chains_1,
            default=available_chains_1,
            key="alignment_chains_1",
        )

    with col2:
        st.markdown("**Chains to Consider in Structure to Align**")
        selected_chains_2 = st.multiselect(
            "Chains to Align",
            available_chains_2,
            default=available_chains_2,
            key="alignment_chains_2",
        )

    return selected_chains_1, selected_chains_2
