import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def plot_real_and_hypothetical(
    coords,
    real_indices,
    option1_vertices=None,
    option2_vertices=None,
    nearest_indices_option1=None,
    nearest_indices_option2=None,
    targetid=None,
):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Real galaxies (red)
    if real_indices is not None and len(real_indices) > 0:
        real_coords = coords[real_indices]
        ax.scatter(
            real_coords[:, 0],
            real_coords[:, 1],
            real_coords[:, 2],
            s=40,
            c="red",
            label="Real galaxies",
        )

    # Option 1 hypothetical vertices (blue)
    if option1_vertices is not None:
        ax.scatter(
            option1_vertices[:, 0],
            option1_vertices[:, 1],
            option1_vertices[:, 2],
            s=60,
            c="blue",
            marker="^",
            label="Hypothetical vertices (Option 1)",
        )

    # Option 2 hypothetical vertices (lighter blue)
    if option2_vertices is not None:
        ax.scatter(
            option2_vertices[:, 0],
            option2_vertices[:, 1],
            option2_vertices[:, 2],
            s=60,
            c="#66aaff",
            marker="s",
            label="Hypothetical vertices (Option 2)",
        )

    # Nearest real galaxies to Option 1 vertices (green)
    if nearest_indices_option1 is not None and len(nearest_indices_option1) > 0:
        ncoords = coords[nearest_indices_option1]
        ax.scatter(
            ncoords[:, 0],
            ncoords[:, 1],
            ncoords[:, 2],
            s=50,
            c="green",
            marker="o",
            label="Nearest to Opt1 vertices",
        )

    # Nearest real galaxies to Option 2 vertices (dark green)
    if nearest_indices_option2 is not None and len(nearest_indices_option2) > 0:
        ncoords2 = coords[nearest_indices_option2]
        ax.scatter(
            ncoords2[:, 0],
            ncoords2[:, 1],
            ncoords2[:, 2],
            s=50,
            c="#008800",
            marker="o",
            label="Nearest to Opt2 vertices",
        )

    ax.set_title("Real and Hypothetical Cube Vertices")
    ax.set_xlabel("X (Mpc)")
    ax.set_ylabel("Y (Mpc)")
    ax.set_zlabel("Z (Mpc)")
    ax.legend(loc="best")

    plt.tight_layout()
    plt.show()
