"""Definition of SPT information."""
from typing import List

from porous_media.visualization.pyvista_visualization import Scalar

scalars_spt: List[Scalar] = [
    Scalar(
        sid="rr_necrosis", title="Necrosis (0: alive, 1: death)", colormap="binary"
    ),
    # Scalar(sid="rr_(atp)", title="ATP (mM)", colormap="RdBu"),
    # Scalar(sid="rr_(glc)", title="Glucose (mM)", colormap="RdBu"),
    # Scalar(sid="rr_(lac)", title="Lactate (mM)", colormap="RdBu"),
    # Scalar(sid="rr_(o2)", title="Oxygen (mM)", colormap="RdBu"),
    # Scalar(sid="rr_(pyr)", title="Pyruvate (mM)", colormap="RdBu"),
    # Scalar(sid="rr_(adp)", title="ADP (mM)", colormap="RdBu"),
    # Scalar(sid="rr_(nadh)", title="NADH (mM)", colormap="RdBu"),
    # Scalar(sid="rr_(nad)", title="NAD (mM)", colormap="RdBu"),
    # Scalar(sid="rr_(ros)", title="ROS (mM)", colormap="RdBu"),
    # Scalar(sid="rr_(alt)", title="ALT (mM)", colormap="RdBu"),
    # Scalar(sid="rr_(ast)", title="AST (mM)", colormap="RdBu"),
    # Scalar(sid="rr_Vext", title="Volume extern (l)", colormap="RdBu"),
    # Scalar(sid="rr_Vli", title="Volume liver (l)", colormap="RdBu"),
]
