"""Plotly viewer components for 3D structure visualization."""

import matplotlib.colors as mcolors
import numpy as np
import plotly.graph_objects as go


def plot_vectors_plotly(
    df,
    cmap,
    vmin,
    vmax,
    show_heads=True,
    line_width=2,
    head_size=3,
    width=1200,
    height=800,
    fidelity=1,
):
    df = df.dropna(subset=["Coord1", "Coord2", "Distance"])
    if df.empty:
        return go.Figure()

    df = df.iloc[::fidelity]

    coords1 = np.array(df["Coord1"].tolist(), dtype=float)
    coords2 = np.array(df["Coord2"].tolist(), dtype=float)
    distances = df["Distance"].astype(float).values

    # Normalize for colormap
    normed = mcolors.Normalize(vmin, vmax)(distances)
    colors_rgb = (np.array([cmap(n)[:3] for n in normed]) * 255).astype(int)
    color_strings = [f"rgb({r},{g},{b})" for r, g, b in colors_rgb]

    x = np.empty((len(coords1) * 3,))
    y = np.empty_like(x)
    z = np.empty_like(x)
    line_colors = []

    for i, (start, end, color) in enumerate(zip(coords1, coords2, color_strings)):
        idx = i * 3
        x[idx : idx + 3] = [start[0], end[0], np.nan]
        y[idx : idx + 3] = [start[1], end[1], np.nan]
        z[idx : idx + 3] = [start[2], end[2], np.nan]
        line_colors.extend([color] * 3)

    fig = go.Figure()

    # Add vectors
    fig.add_trace(
        go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="lines",
            line={"color": line_colors, "width": line_width},
            hoverinfo="skip",
        )
    )

    # Add heads
    if show_heads:
        fig.add_trace(
            go.Scatter3d(
                x=coords2[:, 0],
                y=coords2[:, 1],
                z=coords2[:, 2],
                mode="markers",
                marker={
                    "size": head_size,
                    "color": color_strings,
                    "symbol": "circle",
                },
                hoverinfo="skip",
            )
        )

    # Layout: remove all gridlines, axes, and background
    fig.update_layout(
        scene={
            "xaxis": {
                "showbackground": False,
                "showticklabels": False,
                "showgrid": False,
                "zeroline": False,
                "visible": False,
            },
            "yaxis": {
                "showbackground": False,
                "showticklabels": False,
                "showgrid": False,
                "zeroline": False,
                "visible": False,
            },
            "zaxis": {
                "showbackground": False,
                "showticklabels": False,
                "showgrid": False,
                "zeroline": False,
                "visible": False,
            },
            "aspectmode": "data",
        },
        width=width,
        height=height,
        margin={"l": 0, "r": 0, "t": 0, "b": 0},
        paper_bgcolor="white",
    )

    return fig
