#!/usr/bin/env bash
# One-stop check on Linux/macOS: install, test, run the demo.
set -euo pipefail
python -m pip install -e ".[dev]"
python -m pytest
snowtrack demo --output-dir examples/_work/demo
