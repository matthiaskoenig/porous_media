"""Render FEBio input files from templates and explicit simulation contexts."""

from collections.abc import Mapping
from pathlib import Path

from jinja2 import FileSystemLoader, StrictUndefined
from jinja2.sandbox import ImmutableSandboxedEnvironment

from porous_media._paths import filename_component


def create_feb_files(
    feb_template: Path,
    feb_data: Mapping[str, Mapping[str, object]],
    output_dir: Path,
    *,
    overwrite: bool = False,
) -> None:
    """Render XML-escaped contexts with strict, sandboxed Jinja2 evaluation.

    Args:
        feb_template: UTF-8 Jinja2 FEBio template. Includes resolve beside it.
        feb_data: Template context for each portable simulation identifier.
        output_dir: Directory for generated FEBio files.
        overwrite: Explicitly allow replacing existing regular output files.

    Raises:
        ValueError: An identifier is unsafe or an output is a symbolic link.
        FileExistsError: An output exists and overwrite is false.
        jinja2.TemplateError: Template compilation or rendering fails.

    Templates must still be reviewed and contexts must contain trusted data.
    Sandboxing does not limit CPU, memory, or calls on user-supplied objects.
    """
    prefix = filename_component(feb_template.name.split(".")[0])
    outputs: dict[str, Path] = {}
    for sim_id in feb_data:
        filename_component(sim_id)
        path = output_dir / f"{prefix}_{sim_id}.feb"
        if path.is_symlink():
            raise ValueError(f"Refusing to write through a symbolic link: {path}")
        if path.exists() and not overwrite:
            raise FileExistsError(f"FEBio output already exists: {path}")
        outputs[sim_id] = path

    environment = ImmutableSandboxedEnvironment(
        loader=FileSystemLoader(feb_template.parent, encoding="utf-8"),
        undefined=StrictUndefined,
        autoescape=True,
        trim_blocks=True,
        lstrip_blocks=True,
    )
    template = environment.get_template(feb_template.name)
    for sim_id, context in feb_data.items():
        rendered = template.render(context)
        output_dir.mkdir(parents=True, exist_ok=True)
        with outputs[sim_id].open(
            "w" if overwrite else "x", encoding="utf-8"
        ) as output:
            output.write(rendered)
