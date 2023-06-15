from pathlib import Path
from typing import List, Dict
import jinja2
import numpy as np
import pandas as pd

from pm.console import console


def create_feb_files(feb_template: Path, feb_data: Dict[str, Dict], output_dir: Path):
    """Creates feb files from template by rendering the template tags with information.

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
        console.print(f"{out_path}: {context}")
        with open(out_path, "w") as f_feb:
            f_feb.write(feb_str)


def create_feb_files_xlsx(feb_template: Path, sim_xlsx: Path, output_dir: Path, mapping: Dict[str, str]):
    """Create feb files from the XLSX definition.

    No unit conversion is performed! So the data in excel must be in the right units.

    :param feb_template: path to feb template file with tags
    :param feb_data: dictionary of simulation_id:context information for rendering (provides information for the tags)
    :param output_dir: directory where the feb files are generated.
    :param mapping: mapping of tag:column from excel.
    """
    df = pd.read_excel(sim_xlsx)
    console.print(df)

    feb_data = {}
    for k, row in df.iterrows():
        context = {}
        for tag, column in mapping.items():
            context[tag] = row[column]

        sim_id = row["simulation_id"]
        feb_data[sim_id] = context

    create_feb_files(
        feb_template=feb_template,
        feb_data=feb_data,
        output_dir=feb_out_dir,
    )


if __name__ == "__main__":
    from pm import DATA_DIR
    feb_template = DATA_DIR / "febio" / "lobule_BCflux.feb.template"

    # -----------------------------
    # Example1: simple rendering
    # -----------------------------
    console.rule(title="Example1", align="left", style="white")
    feb_out_dir: Path = DATA_DIR / "febio" / "febs_example1"

    feb_data = {
        "sim1": {
            "boundary_pressure": 133,
            "boundary_flux": 1
        },
        "sim2": {
            "boundary_pressure": 266,
            "boundary_flux": 2
        }
    }
    create_feb_files(
        feb_template=feb_template,
        feb_data=feb_data,
        output_dir=feb_out_dir,
    )

    # --------------------------------
    # Example2: use xlsx information
    # --------------------------------
    console.rule(title="Example2", align="left", style="white")
    feb_template = DATA_DIR / "febio" / "lobule_BCflux.feb.template"
    feb_out_dir: Path = DATA_DIR / "febio" / "scanpf_xlsx"
    xlsx_path: Path = DATA_DIR / "febio" / "scan_pressure_flow.xlsx"
    create_feb_files_xlsx(
        feb_template=feb_template,
        sim_xlsx=xlsx_path,
        output_dir=feb_out_dir,
        mapping={
            "boundary_pressure": "pressure_pf",
            "boundary_flux": "substrate_flow"
        }
    )

    # --------------------------------
    # Example3: use numpy vectors
    # --------------------------------
    console.rule(title="Example3", align="left", style="white")
    boundary_pressure = np.linspace(start=1.37*133, stop=5*133, num=10)
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
    feb_out_dir: Path = DATA_DIR / "febio" / "scanpf"
    create_feb_files(
        feb_template=feb_template,
        feb_data=feb_data,
        output_dir=feb_out_dir,
    )
