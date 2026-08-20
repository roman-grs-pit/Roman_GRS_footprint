"""Tests for the rstgrs_footprint package."""

import pytest

from rstgrs_footprint.example import add_one


def test_add_one_numbers_and_string():
    """add_one should increment numeric inputs and append to strings."""
    assert add_one(1) == 2
    assert add_one(1.1) == 2.1
    assert add_one("1") == "11"


def test_add_one_rejects_unsupported_types():
    """add_one should raise TypeError for unsupported input types."""
    with pytest.raises(TypeError):
        add_one([1])
