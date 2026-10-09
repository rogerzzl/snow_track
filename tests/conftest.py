from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


@pytest.fixture(scope="session")
def examples_dir() -> Path:
    return EXAMPLES


@pytest.fixture(scope="session")
def vic_example(examples_dir) -> Path:
    d = examples_dir / "02_vic_classic_changdu"
    if not d.is_dir():
        pytest.skip("bundled VIC example not available")
    return d


@pytest.fixture(scope="session")
def outlets_dir(examples_dir) -> Path:
    d = examples_dir / "03_basin_outlet_metrics" / "data"
    if not d.is_dir():
        pytest.skip("bundled outlet example not available")
    return d
