# Main package exports
from .core.auto_alignment import StructureMapper, EnsembleMapper
from .core.structure_parsing import load_structure, get_coords_from_id
from .core.distance_calc import DistanceMatrix, CompareDistanceMatrix, calc_difference_from_mapper

from .viz.plotting import plot_distance_difference, plot_interactive_contact_map
from .viz.export_visualizations import (
    get_color_discrete,
    get_color_gradient,
    generate_chimera_link_script,
    generate_arrow_dicts,
    generate_multiple_movement_scripts
)

__all__ = [
    # Core
    'StructureMapper',
    'EnsembleMapper', 
    'load_structure',
    'get_coords_from_id',
    'DistanceMatrix',
    'CompareDistanceMatrix',
    'calc_difference_from_mapper',
    
    # Visualization
    'plot_distance_difference',
    'plot_interactive_contact_map',
    'get_color_discrete',
    'get_color_gradient',
    'generate_chimera_link_script',
    'generate_arrow_dicts',
    'generate_multiple_movement_scripts',
]
