import numpy as np

EDGE = 137.5
TOL = 10.0   # Option B tolerance

# Fixed cube axis
BASE_AXIS = np.array([0.3, 0.8, 0.5], dtype=float)
BASE_AXIS /= np.linalg.norm(BASE_AXIS)


def build_fixed_cube(center, coords, tree):
    """
    Build a cube with:
    - corner 0 = center galaxy
    - fixed orientation using BASE_AXIS
    - edge length = 137.5 Mpc
    - returns 8 vertices (corner positions)
    """

    u = BASE_AXIS

    # Build a perpendicular vector v
    if abs(u[0]) < 0.9:
        temp = np.array([1, 0, 0], float)
    else:
        temp = np.array([0, 1, 0], float)

    v = temp - np.dot(temp, u) * u
    v /= np.linalg.norm(v)

    # Third axis w = u × v
    w = np.cross(u, v)
    w /= np.linalg.norm(w)

    # Scale axes to cube edge
    e1 = EDGE * u
    e2 = EDGE * v
    e3 = EDGE * w

    # Build 8 corners
    vertices = []
    for a in [0, 1]:
        for b in [0, 1]:
            for c in [0, 1]:
                vtx = center + a * e1 + b * e2 + c * e3
                vertices.append(vtx)

    return np.array(vertices)
