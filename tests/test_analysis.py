import pytest
import numpy as np
from pathlib import Path
from src import analysis

def test_moving_average():
    data = [1, 2, 3, 4, 5]
    result = analysis.moving_average(data, 3)
    # Expected: [(1+2+3)/3, (2+3+4)/3, (3+4+5)/3] = [2, 3, 4]
    expected = [2.0, 3.0, 4.0]
    assert np.allclose(result, expected)

def test_moving_average_window_1():
    data = [1, 2, 3]
    result = analysis.moving_average(data, 1)
    assert result == data

def test_calculate_encapsulation_efficiency_empty():
    # Test with non-existent path
    eff = analysis.calculate_encapsulation_efficiency(Path("non_existent.dump"), 1, 1, 1.0)
    assert eff == 0.0

# Mock dump file for testing
@pytest.fixture
def mock_dump_file(tmp_path):
    dump_content = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
4
ITEM: BOX BOUNDS pp pp pp
-10 10
-10 10
-10 10
ITEM: ATOMS id mol type x y z
1 1 1 0.0 0.0 0.0
2 1 1 1.0 1.0 1.0
3 2 2 0.5 0.5 0.5
4 2 2 5.0 5.0 5.0
"""
    dump_path = tmp_path / "test.dump"
    with open(dump_path, "w") as f:
        f.write(dump_content)
    return dump_path

def test_calculate_encapsulation_efficiency_basic(mock_dump_file):
    # Polymer atoms: (0,0,0) and (1,1,1) -> COM = (0.5, 0.5, 0.5)
    # Payload 1 (mol 2): COM = (0.5+5.0)/2 = (2.75, 2.75, 2.75)
    # Distance: sqrt(3 * (2.75 - 0.5)^2) = sqrt(3 * 2.25^2) = 2.25 * 1.732 = 3.897
    # rg = 1.0. threshold = 1.5 * rg = 1.5. Distance 3.897 > 1.5 -> efficiency 0%
    eff = analysis.calculate_encapsulation_efficiency(mock_dump_file, 1, 1, 1.0)
    assert eff == 0.0

    # Increase rg to make it encapsulate
    eff = analysis.calculate_encapsulation_efficiency(mock_dump_file, 1, 1, 5.0)
    assert eff == 100.0
