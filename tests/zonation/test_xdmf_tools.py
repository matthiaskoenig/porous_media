"""Test xdmf functionality."""

import pytest

from porous_media import RESOURCES_DIR
from porous_media.data.xdmf_tools import XDMFInfo, vtks_to_xdmf


def test_vtk_single_to_xdmf(tmp_path):
    """Test the vtk time course parsing."""
    vtk_dir = RESOURCES_DIR / "vtk" / "vtk_single"
    xdmf_path = tmp_path / "vtk_single.xdmf"
    vtks_to_xdmf(vtk_dir, xdmf_path=xdmf_path, overwrite=True)
    assert xdmf_path.exists()

    xdmf_info: XDMFInfo = XDMFInfo.from_path(xdmf_path)
    assert xdmf_info.num_steps == 1
    assert xdmf_info.tstart == pytest.approx(220.0)
    assert xdmf_info.tend == pytest.approx(220.0)

    for key in [
        "displacement",
        "effective_fluid_pressure_TPM",
    ]:
        assert key in xdmf_info.point_data

    for key in [
        "stress",
        "fluid_flux_TPM",
        "solid_stress_TPM",
        "S_ext",
        "P_ext",
    ]:
        assert key in xdmf_info.cell_data


def test_vtk_timecourse_to_xdmf(tmp_path):
    """Test the vtk time course parsing."""
    vtk_dir = RESOURCES_DIR / "vtk" / "vtk_timecourse"
    xdmf_path = tmp_path / "vtk_timecourse.xdmf"
    vtks_to_xdmf(vtk_dir, xdmf_path=xdmf_path, overwrite=True)
    assert xdmf_path.exists()

    xdmf_info: XDMFInfo = XDMFInfo.from_path(xdmf_path)
    assert xdmf_info.num_steps == 3
    assert xdmf_info.tstart == pytest.approx(0.0)
    assert xdmf_info.tend == pytest.approx(220.0)
