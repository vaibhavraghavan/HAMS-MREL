"""
Compute generalized mode shape function for a HAMS hull mesh.

Calculates the normal = nx*X + ny*Y + nz*Z at each panel centroid in global
coordinates, using the Local Coordinate System (LCS) from the HAMS ControlFile.

Usage:
    python panel_centroids.py <mesh_file> <control_file> <dof> <nx> <ny> <nz>

Example:
    python panel_centroids.py HullMesh_1.pnl ControlFile.in 7 0.0 0.0 1.0
"""

import sys
import os
import numpy as np


def read_pnl_mesh(filepath):
    """
    Read a HAMS .pnl hull mesh file and extract node coordinates and panel connectivity.

    Parameters
    ----------
    filepath : str
        Path to the .pnl mesh file.

    Returns
    -------
    nodes : dict
        Dictionary mapping node ID (int) -> np.array([x, y, z])
    panels : list of dict
        Each entry has 'panel_id' (int), 'num_vertices' (int), 'vertex_ids' (list of int)
    n_panels : int
        Number of panels
    n_nodes : int
        Number of nodes
    """
    nodes = {}
    panels = []
    n_panels = 0
    n_nodes = 0

    with open(filepath, 'r') as f:
        lines = f.readlines()

    section = None

    for line in lines:
        stripped = line.strip()

        if not stripped or stripped.startswith('-'):
            continue

        if '# Number of Panels' in stripped:
            section = 'header'
            continue
        elif '# Start Definition of Node Coordinates' in stripped:
            section = 'nodes'
            continue
        elif '# End Definition of Node Coordinates' in stripped:
            section = None
            continue
        elif '# Start Definition of Node Relations' in stripped:
            section = 'panels'
            continue
        elif '# End Definition of Node Relations' in stripped:
            section = None
            continue

        if stripped.startswith('#'):
            continue

        tokens = stripped.split()

        if section == 'header':
            n_panels = int(tokens[0])
            n_nodes = int(tokens[1])
            section = None

        elif section == 'nodes':
            node_id = int(tokens[0])
            x, y, z = float(tokens[1]), float(tokens[2]), float(tokens[3])
            nodes[node_id] = np.array([x, y, z])

        elif section == 'panels':
            panel_id = int(tokens[0])
            num_vertices = int(tokens[1])
            vertex_ids = [int(tokens[2 + k]) for k in range(num_vertices)]
            panels.append({
                'panel_id': panel_id,
                'num_vertices': num_vertices,
                'vertex_ids': vertex_ids,
            })

    return nodes, panels, n_panels, n_nodes


def extract_body_number(mesh_filepath):
    """
    Extract the body number from the mesh filename.
    E.g., 'HullMesh_1.pnl' -> 1, 'HullMesh_2.pnl' -> 2
    """
    basename = os.path.basename(mesh_filepath)
    name_no_ext = os.path.splitext(basename)[0]  # 'HullMesh_1'
    body_number = int(name_no_ext.split('_')[-1])
    return body_number


def read_lcs_from_controlfile(control_filepath, body_number):
    """
    Read the Local Coordinate System (LCS) for a given body from the HAMS ControlFile.

    Parameters
    ----------
    control_filepath : str
        Path to ControlFile.in
    body_number : int
        Body number (1-indexed).

    Returns
    -------
    x0, y0, z0 : float
        Translation offsets (local origin in global frame).
    rot_deg : float
        Rotation angle about the z-axis in degrees.
    """
    lcs_key = f'LCS_{body_number}'

    with open(control_filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith(lcs_key):
                tokens = stripped.split()
                # Format: LCS_N  x0  y0  z0  rot_angle
                x0 = float(tokens[1])
                y0 = float(tokens[2])
                z0 = float(tokens[3])
                rot_deg = float(tokens[4])
                return x0, y0, z0, rot_deg

    raise ValueError(f"LCS for body {body_number} ('{lcs_key}') not found in {control_filepath}")


def local_to_global(coords_local, x0, y0, z0, rot_deg):
    """
    Transform coordinates from local body frame to global frame.

    Applies rotation about the z-axis followed by translation.

    Parameters
    ----------
    coords_local : np.ndarray, shape (N, 3)
        Local coordinates.
    x0, y0, z0 : float
        Translation offsets.
    rot_deg : float
        Rotation angle about z-axis in degrees.

    Returns
    -------
    coords_global : np.ndarray, shape (N, 3)
    """
    rot_rad = np.radians(rot_deg)
    cos_r = np.cos(rot_rad)
    sin_r = np.sin(rot_rad)

    # Rotation matrix about z-axis
    R = np.array([
        [cos_r, -sin_r, 0.0],
        [sin_r,  cos_r, 0.0],
        [0.0,    0.0,   1.0]
    ])

    # Rotate then translate
    coords_global = (R @ coords_local.T).T + np.array([x0, y0, z0])

    return coords_global


def compute_centroids(nodes, panels):
    """
    Compute the centroid of each panel as the arithmetic mean of its vertex coordinates.

    Returns
    -------
    panel_ids : np.ndarray, shape (n_panels,)
    centroids : np.ndarray, shape (n_panels, 3)
    """
    n = len(panels)
    panel_ids = np.zeros(n, dtype=int)
    centroids = np.zeros((n, 3))

    for i, panel in enumerate(panels):
        panel_ids[i] = panel['panel_id']
        coords = np.array([nodes[vid] for vid in panel['vertex_ids']])
        centroids[i] = coords.mean(axis=0)

    return panel_ids, centroids


def compute_generalized_normal(centroids_global, nx, ny, nz):
    """
    Compute the generalized mode shape function at each panel centroid.

    normal_j = nx * X + ny * Y + nz * Z

    Parameters
    ----------
    centroids_global : np.ndarray, shape (N, 3)
        Global centroid coordinates.
    nx, ny, nz : float
        Mode shape direction coefficients.

    Returns
    -------
    normals : np.ndarray, shape (N,)
    """
    return nx * centroids_global[:, 0] + ny * centroids_global[:, 1] + nz * centroids_global[:, 2]


def main():

    mesh_file = 'HullMesh_1.pnl'
    control_file = 'ControlFile.in'
    dof = 1

    # --- Extract body number from mesh filename ---
    body_number = extract_body_number(mesh_file)
    print(f"Mesh file:    {mesh_file}")
    print(f"Control file: {control_file}")
    print(f"Body number:  {body_number}")
    print(f"DOF:          {dof}")

    # --- Read mesh ---
    nodes, panels, n_panels, n_nodes = read_pnl_mesh(mesh_file)
    print(f"\nMesh: {n_panels} panels, {n_nodes} nodes")

    # --- Compute local centroids ---
    panel_ids, centroids_local = compute_centroids(nodes, panels)

    # --- Read LCS and transform to global ---
    x0, y0, z0, rot_deg = read_lcs_from_controlfile(control_file, body_number)
    print(f"LCS_{body_number}: x0={x0}, y0={y0}, z0={z0}, rot={rot_deg} deg")

    centroids_global = local_to_global(centroids_local, x0, y0, z0, rot_deg)

    nx = 1
    ny = 2
    nz = 3

    # --- Compute generalized normal ---
    normals = compute_generalized_normal(centroids_global, nx, ny, nz)
    print(f"nx, ny, nz:   {nx}, {ny}, {nz}")

    # --- Write output ---
    output_file = f"gen_mod_{body_number}_{dof}.txt"
    with open(output_file, 'w') as f:
        for pid, normal in zip(panel_ids, normals):
            f.write(f"  {pid:>6}    {normal:>18.6E}\n")

    print(f"\nOutput written to: {output_file}")
    print(f"  {len(panel_ids)} panels written")
    print(f"\nFirst 10 entries:")
    print(f"{'Panel':>6}  {'Normal':>18}")
    print("-" * 28)
    for pid, normal in zip(panel_ids[:10], normals[:10]):
        print(f"{pid:>6}  {normal:>18.6E}")


if __name__ == '__main__':
    main()
