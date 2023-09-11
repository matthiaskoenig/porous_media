"""Create febio files.

dependencies:
    numpy
    jinja2
"""


from pathlib import Path
from typing import Dict, List

import jinja2
import numpy as np


# from porous_media.console import console


def create_feb_files(
    feb_template: Path, feb_data: Dict[str, Dict], output_dir: Path
) -> None:
    """Create feb files from template by rendering the template tags with information.

    Tags are of the form: {{ boundary_pressure }}, {{ boundary_flux }}
    Uses jinja2 language for context rendering.

    :param feb_template: path to feb template file with tags
    :param feb_data: dictionary of simulation_id:context information for rendering (provides information for the tags)
    :param output_dir: directory where the feb files are generated.
    """
    # create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # template environment
    template_dir = feb_template.parent
    template_file = feb_template.name
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(template_dir),
        extensions=[],
        trim_blocks=True,
        lstrip_blocks=True,
    )
    template = env.get_template(template_file)

    # render all contexts:
    prefix = template_file.split(".")[0]
    for sim_id, context in feb_data.items():
        feb_str = str(template.render(context))
        out_path = output_dir / f"{prefix}_{sim_id}.feb"
        # console.print(f"{out_path}: {context}")
        print(f"{out_path}: {context}")
        with open(out_path, "w") as f_feb:
            f_feb.write(feb_str)


if __name__ == "__main__":
    from porous_media import DATA_DIR

    feb_template = DATA_DIR / "febio" / "lobule_BCflux.feb.template"

    # -----------------------------
    # Example1: simple rendering
    # -----------------------------
    # console.rule(title="Example1", align="left", style="white")
    print("-" * 80)
    feb_out_dir: Path = DATA_DIR / "febio" / "febs_example1"

    feb_data = {
        "sim1": {"boundary_pressure": 133, "boundary_flux": 1},
        "sim2": {"boundary_pressure": 266, "boundary_flux": 2},
    }
    create_feb_files(
        feb_template=feb_template,
        feb_data=feb_data,
        output_dir=feb_out_dir,
    )

    # --------------------------------
    # Example2: use numpy vectors
    # --------------------------------
    # console.rule(title="Example2", align="left", style="white")
    print("-" * 80)
    boundary_pressure = np.linspace(start=1.37 * 133, stop=5 * 133, num=10)
    boundary_flow = np.linspace(start=0.00008441679, stop=0.00000422083951, num=10)

    feb_data = {}
    counter = 1
    for pressure in boundary_pressure:
        for flow in boundary_flow:
            feb_data[f"sim{counter:>03}"] = {
                "boundary_pressure": pressure,
                "boundary_flow": flow,
            }
            counter += 1

    feb_template = DATA_DIR / "febio" / "lobule_BCflux.feb.template"
    feb_out_dir = DATA_DIR / "febio" / "scanpf"
    create_feb_files(
        feb_template=feb_template,
        feb_data=feb_data,
        output_dir=feb_out_dir,
    )
