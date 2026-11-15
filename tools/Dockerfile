# The tools' environment (tools/*.py): the Python packages at the versions of the published results, and the
# repository's NumPy engine (qsim) on the path. Build from the repository root:
#   docker build -t qcs-tools -f tools/Dockerfile .
#   docker run --rm -v "$PWD":/repo -w /repo qcs-tools python tools/<script>.py ...
FROM python:3.12-slim@sha256:78387bc3881b8273120a12ebe6c1ab22b018ccc2c9adf565ae1ac9b536e184ea
# Stim and PyMatching have no wheels for every platform (aarch64 here): pip compiles them
RUN apt-get update && apt-get install -y --no-install-recommends g++ cmake make git ca-certificates && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==2.5.4 scipy==1.18.1 sympy==1.14.0 networkx==3.7 matplotlib==3.11.2 \
        stim==1.16.0 PyMatching==2.4.0 qiskit==2.6.1 qiskit-aer==0.18.2 perceval-quandela==1.3.1 exqalibur==1.4.2 pymupdf==2.0
COPY qsim/qsim /opt/qsim/qsim
ENV PYTHONPATH=/opt/qsim PYTHONDONTWRITEBYTECODE=1
