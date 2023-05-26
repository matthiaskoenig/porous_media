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


def read_simulation_xlsx(base_dir: Path) -> pd.DataFrame:
    """Read the simulation description."""

    df = pd.read_excel(base_dir / "simulation.xlsx", sheet_name=0, comment="#")
    # drop empty lines
    df.dropna(how="all", inplace=True)

    validate_simulation_df(df)
    return df


def validate_simulation_df(df):
    """Validation of simulation DataFrame."""

    columns = set(df.columns)
    # assert simulation id

    # validate lowercase
    for col in columns:
        if col.lower() != col:
            logger.error(f"Columns must be lower case: '{col}'")

    # no whitespaces or special characters
    for col in columns:
        if not re.match("^[a-zA-Z0-9_]*$", col):
            msg = (
                f"No special characters or whitespace allowed in column names: '{col}'. "
                f"Allowed characters are 'a-zA-Z0-9_'."
            )
            logger.error(msg)

    # check naming of unit columns
    suffixes = ["_unit"]
    for suffix in suffixes:
        required_column_names = [
            col[: -len(suffix)]
            for col in columns
            if col.endswith(suffix)
        ]
        for rcol in required_column_names:
            if rcol not in columns:
                msg = (
                    f"If a column was named with the following pattern: "
                    f"<*{suffix}>. Then the column name <*> is required. "
                    f"Column '{rcol + suffix}' was defined but "
                    f"not '{rcol}'."
                )
                logger.error(msg)

    # TODO: check columns order, i.e. unit columns after corresponding columns;
    # for col in columns

    # TODO: check that units can be parsed
    # for col in columns


if __name__ == "__main__":
    EXAMPLE_SIMULATION_XLSX
    base_dir = Path("/data/qualiperf/P7-Perf/spt_results")
    df = read_simulation_xlsx(base_dir)
    console.print(df.to_string())
