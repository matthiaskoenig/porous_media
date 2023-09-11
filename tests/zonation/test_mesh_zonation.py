"""Test that all resources for the app are available."""

from protein_distribution import (
    DATA_MERGED_XLSX,
    DATA_RAW_XLSX,
    DATA_XLSX,
    UNIPROT_PATH,
)


def test_app1() -> None:
    """Test that all resources can be loaded for the app."""
    assert UNIPROT_PATH.exists()


def test_app2() -> None:
    """Test that all resources can be loaded for the app."""
    assert DATA_RAW_XLSX.exists()


def test_app3() -> None:
    """Test that all resources can be loaded for the app."""
    assert DATA_MERGED_XLSX.exists()


def test_app4() -> None:
    """Test that all resources can be loaded for the app."""
    assert DATA_XLSX.exists()
