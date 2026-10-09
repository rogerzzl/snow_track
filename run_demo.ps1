# Windows PowerShell one-stop check: install, test, run the demo.
python -m pip install -e ".[dev]"
python -m pytest
snowtrack demo --output-dir examples/_work/demo
