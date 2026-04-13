"""
reef_slice_analysis.py
======================
Compute per-slice volume and XZ-projected area of an STL model.

Approach
--------
Uses trimesh for robust STL loading and slicing:
  - trimesh.intersections.slice_mesh_plane() clips the mesh to each slab
  - trimesh caps the open slice into a watertight solid automatically
  - .volume on a watertight trimesh gives the exact volume
  - XZ projected area uses Shapely union of projected triangles

Dependencies
------------
    pip install trimesh numpy shapely matplotlib

    Optionally, for better boolean operations:
    pip install manifold3d
"""

import numpy as np
import trimesh
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.ops import unary_union
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def analyse_reef_slices(
    stl_path: str,
    n: int,
    plot: bool = False,
) -> List[Dict]:
    """
    Slice an STL mesh into n horizontal layers and return per-slice
    volume and XZ-projected area.

    Parameters
    ----------
    stl_path : str   Path to binary or ASCII STL.
    n        : int   Number of equal-thickness horizontal slices.
    plot     : bool  Show a bar-chart summary.

    Returns
    -------
    List of dicts (one per slice, bottom to top):
        slice_index       int
        z_min, z_max      float  slice boundaries (mesh units)
        volume            float  units^3
        projected_area_xz float  units^2  (projection onto XZ plane, i.e. along Y axis)
    """
    mesh = trimesh.load(stl_path, force='mesh')

    if not isinstance(mesh, trimesh.Trimesh):
        raise ValueError("STL did not load as a single mesh. "
                         "If you have multiple bodies, export them merged from Onshape.")

    # Fix any minor issues (merge duplicate vertices, fix winding, etc.)
    mesh.merge_vertices()
    mesh.fix_normals()

    z_min = float(mesh.bounds[0][2])
    z_max = float(mesh.bounds[1][2])
    slice_edges = np.linspace(z_min, z_max, n + 1)

    results = []
    for i in range(n):
        z_lo = slice_edges[i]
        z_hi = slice_edges[i + 1]

        vol, area = _compute_slice(mesh, z_lo, z_hi)

        results.append(dict(
            slice_index=i,
            z_min=z_lo,
            z_max=z_hi,
            volume=vol,
            projected_area_xz=area,
        ))

    if plot:
        _plot_results(results)

    return results


# ---------------------------------------------------------------------------
# Per-slice computation
# ---------------------------------------------------------------------------

def _compute_slice(
    mesh: trimesh.Trimesh,
    z_lo: float,
    z_hi: float,
) -> Tuple[float, float]:
    """
    Clip mesh to slab [z_lo, z_hi], compute volume and XZ projected area.
    """
    # Clip to z >= z_lo
    sliced = trimesh.intersections.slice_mesh_plane(
        mesh,
        plane_normal=[0, 0, 1],
        plane_origin=[0, 0, z_lo],
        cap=True,   # adds a flat cap to close the cut face
    )

    if sliced is None or sliced.is_empty:
        return 0.0, 0.0

    # Clip to z <= z_hi  (flip normal direction)
    sliced = trimesh.intersections.slice_mesh_plane(
        sliced,
        plane_normal=[0, 0, -1],
        plane_origin=[0, 0, z_hi],
        cap=True,
    )

    if sliced is None or sliced.is_empty:
        return 0.0, 0.0

    # Volume — trimesh gives this directly for a watertight mesh
    # cap=True above should make it watertight; check and warn if not
    if not sliced.is_watertight:
        # Try to repair
        sliced.fill_holes()
        sliced.fix_normals()

    vol = abs(float(sliced.volume))

    # XZ projected area — project all triangles onto XZ plane, take union
    tris = sliced.triangles   # (M, 3, 3)
    area = _projected_area_xz(tris)

    return vol, area


# ---------------------------------------------------------------------------
# XZ projected area
# ---------------------------------------------------------------------------

def _projected_area_xz(tris: np.ndarray) -> float:
    """
    Project triangles onto the XZ plane (drop Y) and return union area.
    This is the silhouette seen by an observer looking along the +Y axis.
    """
    if len(tris) == 0:
        return 0.0

    polys = []
    for tri in tris:
        pts = [(float(v[0]), float(v[2])) for v in tri]
        try:
            p = ShapelyPolygon(pts)
            if p.is_valid and not p.is_empty and p.area > 0:
                polys.append(p)
        except Exception:
            pass

    if not polys:
        return 0.0

    return float(unary_union(polys).area)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def _plot_results(results: List[Dict]) -> None:
    indices = [r['slice_index'] for r in results]
    z_mids  = [(r['z_min'] + r['z_max']) / 2 for r in results]
    vols    = [r['volume'] for r in results]
    areas   = [r['projected_area_xz'] for r in results]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].barh(indices, vols, color='steelblue')
    axes[0].set_yticks(indices)
    axes[0].set_yticklabels([f"{z:.3f}" for z in z_mids])
    axes[0].set_xlabel("Volume (units^3)")
    axes[0].set_ylabel("Slice mid-Z")
    axes[0].set_title("Volume per slice")

    axes[1].barh(indices, areas, color='coral')
    axes[1].set_yticks(indices)
    axes[1].set_yticklabels([f"{z:.3f}" for z in z_mids])
    axes[1].set_xlabel("Projected Area - XZ (units^2)")
    axes[1].set_title("XZ Projected Area per slice")

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("Usage: python reef_slice_analysis.py <file.stl> <n_slices> [--plot]")
        sys.exit(1)

    stl_file = sys.argv[1]
    n_slices = int(sys.argv[2])
    do_plot  = "--plot" in sys.argv

    data = analyse_reef_slices(stl_file, n=n_slices, plot=do_plot)

    print(f"\n{'Slice':>5}  {'z_min':>10}  {'z_max':>10}  "
          f"{'Volume':>14}  {'ProjArea_XZ':>14}")
    print("-" * 62)
    for r in data:
        print(f"{r['slice_index']:>5d}  {r['z_min']:>10.4f}  {r['z_max']:>10.4f}  "
              f"{r['volume']:>14.6f}  {r['projected_area_xz']:>14.6f}")

    total_vol  = sum(r['volume']            for r in data)
    total_area = sum(r['projected_area_xz'] for r in data)
    print("-" * 62)
    print(f"{'TOTAL':>5}  {'':>10}  {'':>10}  {total_vol:>14.6f}  {total_area:>14.6f}")
    print()
    expected = data[0]['volume'] if data else 0
    print(f"All slices equal volume? "
          f"{'YES' if all(abs(r['volume'] - expected) < 1e-6 for r in data) else 'NO — check output above'}")