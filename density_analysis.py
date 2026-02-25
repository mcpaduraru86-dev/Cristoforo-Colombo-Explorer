import numpy as np

def analyze_corner(vtx, coords, tree, targetids):
    """
    For a single cube corner:
    - find nearest galaxy
    - compute distance
    - compute density within 10 Mpc
    """

    dist, idx = tree.query(vtx)
    nearest_id = targetids[idx]

    # Count neighbors within 10 Mpc
    neigh = tree.query_ball_point(vtx, 10.0)
    density = len(neigh) / (4/3 * np.pi * 10.0**3)

    if dist < 10.0:
        label = "GALAXY CORNER"
    elif density > 0.01:
        label = "DENSE REGION"
    else:
        label = "VOID CORNER"

    return {
        "nearest_id": nearest_id,
        "nearest_dist": dist,
        "density": density,
        "label": label
    }
