# Convenience targets (GNU make). On Windows use run_demo.ps1 or the snowtrack CLI directly.
.PHONY: install test demo reproduce clean

install:
	python -m pip install -e ".[dev]"

test:
	python -m pytest

demo:
	snowtrack demo --output-dir examples/_work/demo

reproduce:
	python examples/02_vic_classic_changdu/run_example.py
	python examples/03_basin_outlet_metrics/basin_outlet_metrics.py

clean:
	rm -rf examples/_work examples/*/output .pytest_cache
