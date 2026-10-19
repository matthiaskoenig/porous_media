# Visualization

See [visualization and video](../visualization.md) for rendering requirements,
field selection, camera settings, fixed color limits, and output naming.

## Mesh rendering

::: porous_media.visualization.pyvista_visualization
    options:
      members:
        - DataRangeType
        - DataLayer
        - VisualizationSettings
        - visualize_data_layers
        - visualize_datalayers_timecourse
        - visualize_interactive
        - xdmf_to_mesh
        - create_combined_images

## Image processing

::: porous_media.visualization.image_manipulation
    options:
      members:
        - merge_images

::: porous_media.visualization.image_annotation
    options:
      members:
        - annotate_image_text

## Video

These helpers require FFmpeg on `PATH`.

::: porous_media.visualization.video
    options:
      members:
        - create_video
        - create_gif_from_video
