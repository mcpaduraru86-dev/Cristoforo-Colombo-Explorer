import numpy as np
from scipy.spatial import cKDTree

EDGE = 137.5
FACE_DIAG = EDGE * np.sqrt(2)
SPACE_DIAG = EDGE * np.sqrt(3)

EDGE_TOL = 0.03 * EDGE
FACE_TOL = 0.03 * FACE_DIAG
SPACE_TOL = 0.03 * SPACE_DIAG

ANGLE_TOL_DEG = 5

MAX_NEIGHBORS = 40        # hard limit to avoid explosion
TIME_LIMIT = 2.0          # seconds per galaxy


def neighbors_at_distance(idx, coords, tree, target=EDGE, tol=EDGE_TOL):
    center = coords[idx]
    rmax = target + tol
    cand = tree.query_ball_point(center, rmax)
    cand = [j for j in cand if j != idx]
    if not cand:
        return np.array([], dtype=int)
    d = np.linalg.norm(coords[cand] - center, axis=1)
    mask = np.abs(d - target) < tol
    neigh = np.array(cand)[mask]

    # Limit neighbors to avoid combinatorial explosion
    if len(neigh) > MAX_NEIGHBORS:
        neigh = neigh[:MAX_NEIGHBORS]

    return neigh


def is_orthogonal(i, j, k, coords):
    v1 = coords[j] - coords[i]
    v2 = coords[k] - coords[i]
    dot = np.dot(v1, v2)
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denom == 0:
        return False
    ang = np.degrees(np.arccos(dot / denom))
    return abs(ang - 90) < ANGLE_TOL_DEG


def find_cubes_for_index(i, coords, tree):
    import time
    start = time.time()

    cubes = []
    neigh = neighbors_at_distance(i, coords, tree)

    if len(neigh) < 3:
        return cubes

    for j in neigh:
        if time.time() - start > TIME_LIMIT:
            return cubes

        for k in neigh:
            if j >= k:
                continue
            if not is_orthogonal(i, j, k, coords):
                continue

            v1 = coords[j] - coords[i]
            v2 = coords[k] - coords[i]
            pred_face = coords[i] + v1 + v2

            idx4 = tree.query(pred_face, distance_upper_bound=EDGE_TOL)[1]
            if idx4 >= len(coords):
                continue

            d_jk = np.linalg.norm(coords[j] - coords[k])
            if abs(d_jk - FACE_DIAG) > FACE_TOL:
                continue

            for m in neigh:
                if m in (j, k):
                    continue
                if not (is_orthogonal(i, j, m, coords) and is_orthogonal(i, k, m, coords)):
                    continue

                v3 = coords[m] - coords[i]

                pred2 = coords[i] + v1 + v3
                pred3 = coords[i] + v2 + v3
                pred4 = coords[i] + v1 + v2 + v3

                idx5 = tree.query(pred2, distance_upper_bound=EDGE_TOL)[1]
                idx6 = tree.query(pred3, distance_upper_bound=EDGE_TOL)[1]
                idx7 = tree.query(pred4, distance_upper_bound=EDGE_TOL)[1]

                if max(idx5, idx6, idx7) >= len(coords):
                    continue

                cube = sorted(set([i, j, k, idx4, m, idx5, idx6, idx7]))
                if len(cube) == 8 and cube not in cubes:
                    cubes.append(cube)

    return cubes
