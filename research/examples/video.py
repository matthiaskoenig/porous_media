"""Archived example; configure research paths before running."""

from porous_media.visualization.video import create_gif_from_video

if __name__ == "__main__":
    from research import RESULTS_DIR

    video_path = (
        RESULTS_DIR / "spt_zonation_patterns" / "simulation_pattern0_100_28800.mp4"
    )
    gif_path = (
        RESULTS_DIR / "spt_zonation_patterns" / "simulation_pattern0_100_28800.gif"
    )
    create_gif_from_video(
        video_path=video_path,
        gif_path=gif_path,
    )
