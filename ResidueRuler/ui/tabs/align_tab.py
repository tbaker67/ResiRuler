"""Alignment tab for US-align structure alignment."""

import io
import os

import numpy as np
import py3Dmol
import streamlit as st
from Bio.PDB import MMCIFIO, MMCIFParser
from src.resiruler.core.auto_alignment import filter_and_write_aligned_maps
from src.resiruler.wrappers.usalign_wrapper import run_usalign_matrix_only

from ui.widgets.utils import (
    create_downloadable_zip,
    full_aligner_ui,
    get_threshold,
    reset_downstream,
    save_temp_file,
    select_chains_for_alignment_ui,
)


def apply_transform(structure, rotation, translation):
    """
    Apply rotation (3x3) and translation (3,) to all atoms in a Bio.PDB Structure,
    including HETATMs, waters, ligands, and all alternate locations.
    """
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # Check for disordered atoms (alternate locations)
                    if atom.is_disordered():
                        for alt_atom in atom.child_dict.values():
                            coord = np.asarray(alt_atom.coord, dtype=float)
                            new_coord = np.dot(rotation, coord) + translation
                            alt_atom.set_coord(new_coord)
                    else:
                        coord = np.asarray(atom.coord, dtype=float)
                        new_coord = np.dot(rotation, coord) + translation
                        atom.set_coord(new_coord)
    return structure


def start_pymol_viewer(reference_cif_str, aligned_cif_str):
    view = py3Dmol.view(width=800, height=600)

    # Add reference structure (model 0) - yellow
    view.addModel(reference_cif_str, "cif")
    view.setStyle({"model": 0}, {"cartoon": {"color": "#115A99"}})

    # Add aligned structure (model 1) - semi-transparent blue
    view.addModel(aligned_cif_str, "cif")
    view.setStyle({"model": 1}, {"cartoon": {"color": "#FFC107"}})

    view.zoomTo()
    return view


def show_align_tab():
    st.header("Run US-align Structure Alignment")

    st.markdown(
        "Upload two structure files (.cif) to run US-align and get the rotation & translation matrices."
    )

    struct1 = st.file_uploader("Upload Reference Structure (Structure 1)", type=["cif"])
    struct2 = st.file_uploader(
        "Upload Structure to Align to Reference (Structure 2)", type=["cif"]
    )

    st.session_state.setdefault("usalign_R", None)
    st.session_state.setdefault("usalign_t", None)
    st.session_state.setdefault("aligned_cif", None)
    st.session_state.setdefault("ref_cif", None)

    if struct1 and struct2:
        reference_structure_path = save_temp_file(struct1)
        structure_to_align_path = save_temp_file(struct2)

        reference_chains, chains_to_align = select_chains_for_alignment_ui(
            reference_structure_path, structure_to_align_path
        )
        reference_chains_str = ",".join(reference_chains)
        chains_to_align_str = ",".join(chains_to_align)

    if st.button("Run US-align"):
        if not all([struct1, struct2]):
            st.error("Please Upload Both Structure Files.")
            return

        try:
            # a new run replaces the reference/aligned CIF strings and the
            # (possibly 4x-duplicated) filtered outputs from a prior run -
            # free those before generating the new ones
            reset_downstream(
                "usalign_R",
                "usalign_t",
                "aligned_cif",
                "ref_cif",
                "filtered_outputs",
            )

            # Align structure 2 to structure 1 (mobile -> ref)
            R, t = run_usalign_matrix_only(
                str(structure_to_align_path),
                str(reference_structure_path),
                reference_chains_str,
                chains_to_align_str,
            )

            st.session_state.usalign_R = R
            st.session_state.usalign_t = t

            # Load original reference structure as string
            ref_str = struct1.getvalue().decode("utf-8")
            st.session_state.ref_cif = ref_str

            # Load structure 2 and apply transform
            parser = MMCIFParser(QUIET=True)
            structure_id = os.path.basename(structure_to_align_path).split(".")[0]
            mobile_structure = parser.get_structure(
                structure_id, str(structure_to_align_path)
            )
            aligned_structure = apply_transform(mobile_structure, R, t)

            # Write aligned structure to string buffer
            io_buffer = io.StringIO()
            io_writer = MMCIFIO()
            io_writer.set_structure(aligned_structure)
            io_writer.save(io_buffer)
            aligned_cif_str = io_buffer.getvalue()
            st.session_state.aligned_cif = aligned_cif_str

            st.success("US-align finished!")

        except Exception as e:
            st.error(f"US-align failed: {e}")

    if (
        st.session_state.get("usalign_R") is not None
        and st.session_state.get("aligned_cif") is not None
    ):
        # Show viewer with both reference and aligned structures
        st.subheader("3D Viewer: Reference (Blue) vs Aligned (Yellow)")

        viewer = start_pymol_viewer(
            st.session_state.ref_cif, st.session_state.aligned_cif
        )
        # Embed the viewer's javascript/HTML into Streamlit
        html = viewer._make_html()
        st.components.v1.html(html, height=650)

        # Download aligned CIF
        st.subheader("Download Options")
        cif_filename = st.text_input(
            "CIF Output Filename", value="aligned_structure.cif"
        )

        st.download_button(
            label="Download Aligned CIF",
            data=st.session_state.aligned_cif,
            file_name=cif_filename,
            mime="chemical/x-mm-cif",
        )
        protein_aligner, nucleotide_aligner, rmsd_scaling_factor = full_aligner_ui(key="align")
        pct_id_threshold = get_threshold(
            "Set a Minimum Percent Identity Threshold for Matching Chains Together",
            "95.0",
            "align_pct_id",
        )

        if st.button("Clean Alignment To Show Only Matched Residues/Chains"):
            if not all([struct1, struct2]):
                st.error("Please upload both structure files.")
            else:
                with st.spinner("Filtering matched residues and chains..."):
                    try:
                        # drop any previously-filtered CIF strings before
                        # generating a new set with the current threshold
                        reset_downstream("filtered_outputs")

                        # Use raw file streams
                        ref_cif_stream = io.StringIO(st.session_state.ref_cif)
                        tgt_cif_stream = io.StringIO(st.session_state.aligned_cif)

                        (
                            ref_chain_cif,
                            ref_res_cif,
                            tgt_chain_cif,
                            tgt_res_cif,
                        ) = filter_and_write_aligned_maps(
                            ref_cif_stream,
                            tgt_cif_stream,
                            protein_aligner,
                            nucleotide_aligner,
                            identity_threshold=pct_id_threshold,
                            rmsd_scale_factor=rmsd_scaling_factor
                        )

                        st.session_state["filtered_outputs"] = {
                            "reference_chain_filtered": ref_chain_cif,
                            "reference_residue_filtered": ref_res_cif,
                            "comparison_chain_filtered": tgt_chain_cif,
                            "comparison_residue_filtered": tgt_res_cif,
                        }

                        st.success("Filtered structures generated!")
                        st.subheader("3D Viewer: Filtered Structures")
                        viewer = start_pymol_viewer(
                            st.session_state["filtered_outputs"][
                                "reference_chain_filtered"
                            ],
                            st.session_state["filtered_outputs"][
                                "comparison_chain_filtered"
                            ],
                        )
                        html = viewer._make_html()
                        st.components.v1.html(html, height=650)

                        st.subheader("Download Filtered Structures")
                        files_dict = {
                            f"{label}.cif": cif_str
                            for label, cif_str in st.session_state[
                                "filtered_outputs"
                            ].items()
                        }
                        zip_buffer = create_downloadable_zip(files_dict)
                        st.download_button(
                            label="Download All Filtered Structures (ZIP)",
                            data=zip_buffer,
                            file_name="filtered_structures.zip",
                            mime="application/zip",
                        )

                    except Exception as e:
                        st.error(f"Filtering failed: {e}")
