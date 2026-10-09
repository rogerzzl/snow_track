# Reproducible environment for snowtrack.
#   docker build -t snowtrack .
#   docker run --rm snowtrack snowtrack demo
#   docker run --rm -v "$PWD/mydata:/data" snowtrack snowtrack run --input-dir /data/in --output-dir /data/out
FROM python:3.11-slim
WORKDIR /opt/snowtrack
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
COPY configs ./configs
COPY examples ./examples
COPY tests ./tests
RUN python -m pip install --no-cache-dir -e ".[all]" pytest
CMD ["snowtrack", "demo"]
