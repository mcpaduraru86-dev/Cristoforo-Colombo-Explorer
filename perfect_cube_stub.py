import numpy as np

def edge_lengths(vertices):
    """Return all 12 edge lengths of a cube defined by 8 vertices."""
    edges = []
    for i in range(8):
        for j in range(i + 1, 8):
            d = np.linalg.norm(vertices[i] - vertices[j])
            edges.append(d)
    edges = np.sort(edges)
    return edges[:12]  # shortest 12 distances


def angle_between(v1, v2):
    """Return angle in degrees between two vectors."""
    dot = np.dot(v1, v2)
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denom == 0:
        return None
    cosang = np.clip(dot / denom, -1.0, 1.0)
    return np.degrees(np.arccos(cosang))


def cube_angle_deviations(vertices):
    """Compute all angles between edges meeting at each vertex."""
    angles = []
    for i in range(8):
        dists = np.linalg.norm(vertices - vertices[i], axis=1)
        idxs = np.argsort(dists)[1:4]  # 3 nearest neighbors
        v1 = vertices[idxs[0]] - vertices[i]
        v2 = vertices[idxs[1]] - vertices[i]
        v3 = vertices[idxs[2]] - vertices[i]
        for a, b in [(v1, v2), (v1, v3), (v2, v3)]:
            ang = angle_between(a, b)
            if ang is not None:
                angles.append(ang)
    return np.array(angles)


def compare_to_perfect_cube(vertices):
    """
    Compare a reconstructed cube to a mathematically perfect cube.
    Returns:
    - mean_edge
    - edge_rms_error
    - angle_rms_error
    - cube_quality_score (0–100)
    """

    if vertices is None or len(vertices) != 8:
        return {"error": "Invalid vertex set"}

    edges = edge_lengths(vertices)
    mean_edge = np.mean(edges)
    edge_rms = np.sqrt(np.mean((edges - mean_edge) ** 2))

    angles = cube_angle_deviations(vertices)
    angle_rms = np.sqrt(np.mean((angles - 90.0) ** 2))

    score = 100.0 - (edge_rms * 5 + angle_rms * 1.5)
    score = max(0.0, min(100.0, score))

    return {
        "mean_edge": float(mean_edge),
        "edge_rms_error": float(edge_rms),
        "angle_rms_error": float(angle_rms),
        "cube_quality_score": float(score),
    }
