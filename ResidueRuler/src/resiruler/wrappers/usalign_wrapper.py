import subprocess
import numpy as np
import shutil
import os

def find_usalign_executable():
    '''
    Assuming usalign installed via conda install -c bioconda usalign
    '''
    # Try to find USalign in current PATH
    usalign_path = shutil.which("USalign")
    if usalign_path:
        return usalign_path

    # If not found, try common conda env locations (adjust as needed)
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidate = os.path.join(conda_prefix, "bin", "USalign")
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate

    # Could add other fallback logic here if needed

    raise FileNotFoundError("USalign executable not found in PATH or CONDA_PREFIX/bin.")


def run_usalign_matrix_only(structure_to_align_path, reference_structure_path):
    usalign_binary = find_usalign_executable()

    cmd = [
        usalign_binary,
        structure_to_align_path,
        reference_structure_path,
        "-m", "-",
        "-mm", "1",
        "-ter", "0"

    ]

    # Run US-align
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"USalign failed with error:\n{result.stderr}")

    R, t = parse_matrix_from_stdout(result.stdout)
    return R, t


def parse_matrix_from_stdout(stdout):
    lines = stdout.splitlines()
    matrix_lines = []

    found_header = False
    for line in lines:
        if "rotation matrix" in line.lower():
            found_header = True
            continue
        if found_header:
            if line.strip() and line.strip()[0] in ("0", "1", "2"):
                matrix_lines.append(line.strip())
            if len(matrix_lines) == 3:
                break

    if len(matrix_lines) != 3:
        raise ValueError("Could not parse rotation matrix from US-align output.")

    R = []
    t = []
    for line in matrix_lines:
        parts = line.split()
        t.append(float(parts[1]))
        R.append([float(parts[2]), float(parts[3]), float(parts[4])])

    return np.array(R), np.array(t)