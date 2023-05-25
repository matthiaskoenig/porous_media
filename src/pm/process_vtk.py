"""Scripts for processing simulation results."""
from pathlib import Path
import re

import pandas as pd
from pm.console import console
from pm.log import get_logger

# FIXME: minor changes simulation file
# -> Simulation_id -> simulation_id
# -> mmHG -> mmHg
# -> reorder columns, so that units are corresponding to other columns

# FIXME: IRI model; scripts
# FIXME: rescale;



logger = get_logger("__name__")

# Mesh:

def read_simulation_xlsx(base_dir: Path) -> pd.DataFrame:
    """Read the simulation description."""

    df = pd.read_excel(base_dir / "simulation.xlsx", sheet_name=0, comment="#")
    # drop empty lines
    df.dropna(how="all", inplace=True)
    return df


def validate_simulation_df(df):
    """Validation of simulation DataFrame."""

    columns = set(df.columns)
    # assert simulation id

    # validate lowercase
    for col in columns:
        if col.lower != col:
            logger.error(f"Columns must be lower case: {col}")

    # no whitespaces or special characters
    for col in columns:

        if not re.match("^[a-zA-Z0-9_]*$", col):
            msg = (
                f"No special characters or whitespace allowed in column names: '{col}'"
                f"Allowed characters are 'a-zA-Z0-9_'."
            )
            logger.error(msg)

    # check naming of unit columns
    suffixes = ["_unit"]
    correct_column_name = True
    for suffix in suffixes:
        required_column_names = [
            col[: -len(suffix)]
            for col in column_names
            if column_name.endswith(suffix)
        ]
        for column_name in required_column_names:
            if column_name not in column_names:
                msg = (
                    f"If a column was named with the following pattern: "
                    f"<*{suffix}>. Then the column name <*> is required. "
                    f"In sheet <{sheet_name}> of file <{study_id}.xlsx> a "
                    f"column <{column_name + suffix}> was defined but "
                    f"not <{column_name}>."
                )
                logger.error(msg)
                correct_column_name = False


    # check that units can be parsed


    # ""Validate that no empty lines exist in sheet."""
    if len() < len(df.index):
        msg = (
            f"No empty lines"
            f"empty lines. Remove the line or add a '#' as the first character in the line."
        )
        logger.warning(msg)


if __name__ == "__main__":
    base_dir = Path("/data/qualiperf/P7-Perf/spt_results")
    df = read_simulation_xlsx(base_dir)
    console.print(df.to_string())
