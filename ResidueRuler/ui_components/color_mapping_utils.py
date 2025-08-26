import streamlit as st
import streamlit.components.v1 as components
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap


def get_coloring_values(min=0, max=100):
    """
    UI widget for choosing min/max value.
    """

    col1, col2 = st.columns(2)
    with col1:
        min_val = st.number_input("Min value", value=min)
    with col2:
        max_val = st.number_input("Max value", value=max)

    if min_val >= max_val:
        st.error("Min value must be less than max value")

    return min_val, max_val


def discrete_palette_picker(label="Discrete Color Thresholds", default_thresholds=None, default_colors=None):
    """
    UI to define a list of (threshold, color) mappings for discrete color assignment.
    Returns: list of (threshold: float, color: hex str)
    """
    st.subheader(label)

    num_stops = st.number_input(
        "Number of thresholds", min_value=1, max_value=10,
        value=len(default_thresholds) if default_thresholds else 3,
        step=1, key="num_thresholds"
    )

    thresholds_colors = []

    for i in range(num_stops):
        col1, col2 = st.columns([2, 1])
        with col1:
            default_val = default_thresholds[i] if default_thresholds and i < len(default_thresholds) else 0.0
            threshold = st.number_input(f"Threshold {i+1}", value=default_val, key=f"threshold_{i}")
        with col2:
            default_color = default_colors[i] if default_colors and i < len(default_colors) else "#000000"
            color = st.color_picker(f"Color {i+1}", value=default_color, key=f"color_{i}", label_visibility="collapsed")

        thresholds_colors.append((threshold, color))

    return thresholds_colors

def gradient_palette_picker(label="Gradient Palette", default_colors=None, key="palette_picker"):
    st.subheader(label)

    num_stops = st.number_input("Number of gradient stops", min_value=2, max_value=10,
                                value=len(default_colors) if default_colors else 3, step=1)

    colors = []
    for i in range(num_stops):
        default = default_colors[i % len(default_colors)] if default_colors else "#000000"
        color = st.color_picker(f"Color stop {i + 1}", value=default)
        colors.append(color)


    palette = [hex_color for hex_color in colors]

    return palette

def show_gradient_bar(palette, min_val, max_val):
    """
    Render gradient color bar as HTML
    """
    gradient_colors = ", ".join([hex for hex in palette])
    gradient_css = f"linear-gradient(to right, {gradient_colors})"

    gradient_html = f"""
        <div style="
            width: 100%;
            height: 30px;
            background: {gradient_css};
            border: 1px solid #ddd;
            border-radius: 4px;
            margin-bottom: 8px;
        "></div>
        <div style="display: flex; justify-content: space-between; font-size: 12px;">
            <span>{min_val}</span>
            <span>{max_val}</span>
        </div>
    """
    st.markdown(gradient_html, unsafe_allow_html=True)


def show_discrete_bar(discrete_mapping, min_val, max_val):
    """
    Render discrete color bar from value thresholds
    """
    segments_html = ""
    total_range = max_val - min_val if max_val > min_val else 1
    prev_val = min_val

    widths = []
    for val, _ in discrete_mapping[:-1]:
        width = max(1, (val - prev_val) / total_range * 100)
        widths.append(width)
        prev_val = val

    used_width = sum(widths)
    last_val, last_color = discrete_mapping[-1]
    last_width = max(1, 100 - used_width)

    prev_val = min_val
    for i, (val, color) in enumerate(discrete_mapping[:-1]):
        segments_html += f"""
        <div style="
            width: {widths[i]}%;
            height: 30px;
            background-color: {color};
            display: inline-block;
            float: left;
            margin: 0;
            padding: 0;
        "></div>
        """
        prev_val = val

    segments_html += f"""
    <div style="
        width: {last_width}%;
        height: 30px;
        background-color: {last_color};
        display: inline-block;
        float: left;
        margin: 0;
        padding: 0;
    "></div>
    """

    bar_html = f"""
    <div style="
        border: 1px solid #ddd;
        border-radius: 4px;
        overflow: hidden;
        white-space: nowrap;
        width: 100%;
        height: 30px;
        font-size: 0;
    ">
        {segments_html}
        <div style="clear: both;"></div>
    </div>
    <div style="display: flex; justify-content: space-between; font-size: 12px; margin-top: 4px;">
        <span>{min_val}</span>
        <span>{max_val}</span>
    </div>
    """
    components.html(bar_html, height=50)


def build_gradient_cmap(palette, vmin, vmax):
    """
    Builds a matplot color map based on the provided palette which will evenly space out each color across the range [vmin, vmax] after normalized to [0,1]
    """
    positions = [vmin + i * (vmax - vmin) / (len(palette) - 1) for i in range(len(palette))]
    positions = [(pos - vmin)/(vmax - vmin) for pos in positions]  # matplot requires normalized range to [0,1]

    
    pos_rgb = list(zip(positions, palette))

    cmap = LinearSegmentedColormap.from_list("custom_cmap", pos_rgb)
    return cmap


def sort_discrete_mapping(mapping):
    return sorted(mapping, key=lambda x: x[0])
