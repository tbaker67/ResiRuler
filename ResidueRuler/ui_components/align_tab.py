# align_tab.py
import streamlit as st
import numpy as np
import io
import os
import py3Dmol
from Bio.PDB import MMCIFParser, MMCIFIO
from ui_components.utils import save_temp_file
from usalign_wrapper import run_usalign_matrix_only
from src.resiruler.auto_alignment import filter_and_write_aligned_maps

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

    # Add reference structure (model 0) - blue
    view.addModel(reference_cif_str, 'cif')
    view.setStyle({'model': 0}, {'cartoon': {'color': 'blue'}})

    # Add aligned structure (model 1) - semi-transparent red
    view.addModel(aligned_cif_str, 'cif')
    view.setStyle({'model': 1}, {'cartoon': {'color': 'red'}})

    view.zoomTo()
    return view

def show_align_tab():
    st.header("Run US-align Structure Alignment")

    st.markdown("Upload two structure files (.pdb or .cif) to run US-align and get the rotation & translation matrices.")

    struct1 = st.file_uploader("Upload Reference Structure (Structure 1)", type=["pdb", "cif"])
    struct2 = st.file_uploader("Upload Structure to Align (Structure 2)", type=["pdb", "cif"])

    st.session_state.setdefault("usalign_R", None)
    st.session_state.setdefault("usalign_t", None)
    st.session_state.setdefault("aligned_cif", None)
    st.session_state.setdefault("ref_cif", None)

    if st.button("Run US-align"):
        if not all([struct1, struct2]):
            st.error("Please upload both structure files.")
            return

        path1 = save_temp_file(struct1)
        path2 = save_temp_file(struct2)

        try:
            # Align structure 2 to structure 1 (mobile -> ref)
            R, t = run_usalign_matrix_only(str(path2), str(path1))

            st.session_state.usalign_R = R
            st.session_state.usalign_t = t

            # Load original reference structure as string
            ref_str = struct1.getvalue().decode("utf-8")
            st.session_state.ref_cif = ref_str

            # Load structure 2 and apply transform
            parser = MMCIFParser(QUIET=True)
            structure_id = os.path.basename(path2).split('.')[0]
            mobile_structure = parser.get_structure(structure_id, str(path2))
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

    if st.session_state.get("usalign_R") is not None and st.session_state.get("aligned_cif") is not None:
        # Show viewer with both reference and aligned structures
        st.subheader("3D Viewer: Reference (blue) vs Aligned (red)")

        viewer = start_pymol_viewer(st.session_state.ref_cif, st.session_state.aligned_cif)
        # Embed the viewer's javascript/HTML into Streamlit
        html = viewer._make_html()
        st.components.v1.html(html, height=650)

        # Download aligned CIF
        st.subheader("Download Options")
        cif_filename = st.text_input("CIF Output Filename", value="aligned_structure.cif")

        st.download_button(
            label="Download Aligned CIF",
            data=st.session_state.aligned_cif,
            file_name=cif_filename,
            mime="chemical/x-mm-cif"
        )
    
    if st.button("Clean Alignment To Show Only Matched Residues/Chains"):
        if not all([struct1, struct2]):
            st.error("Please upload both structure files.")
        else:
            with st.spinner("Filtering matched residues and chains..."):
                try:
                    # Use raw file streams
                    ref_cif_stream = io.StringIO(st.session_state.ref_cif)
                    tgt_cif_stream = io.StringIO(st.session_state.aligned_cif)

                    ref_chain_cif, ref_res_cif, tgt_chain_cif, tgt_res_cif = filter_and_write_aligned_maps(
                        ref_cif_stream, tgt_cif_stream, identity_threshold=95.0
                    )

                    st.session_state["filtered_outputs"] = {
                        "ref_chain": ref_chain_cif,
                        "ref_res": ref_res_cif,
                        "tgt_chain": tgt_chain_cif,
                        "tgt_res": tgt_res_cif
                    }

                    st.success("Filtered structures generated!")

                    st.subheader("Download Filtered Structures")

                    for label, cif_str in st.session_state["filtered_outputs"].items():
                        st.download_button(
                            label=f"Download {label.replace('_', ' ').title()} CIF",
                            data=cif_str,
                            file_name=f"{label}.cif",
                            mime="chemical/x-mm-cif"
                        )

                    st.subheader("3D Viewer: Filtered Structures")
                    viewer = py3Dmol.view(width=800, height=600)
                    viewer.addModel(st.session_state["filtered_outputs"]["ref_res"], 'cif')
                    viewer.setStyle({'model': 0}, {'cartoon': {'color': 'blue'}})
                    viewer.addModel(st.session_state["filtered_outputs"]["tgt_res"], 'cif')
                    viewer.setStyle({'model': 1}, {'cartoon': {'color': 'red'}})
                    viewer.zoomTo()

                    html = viewer._make_html()
                    st.components.v1.html(html, height=650)

                except Exception as e:
                    st.error(f"Filtering failed: {e}")


    
        

