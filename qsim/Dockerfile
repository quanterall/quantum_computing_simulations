# qsim: quanterall.com's quantum demo API on a NumPy state-vector engine.
FROM python:3.12-slim AS base
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY qsim qsim

# docker build --target test .   -> the tests run during the build; a failure stops it
FROM base AS test
COPY requirements-test.txt .
RUN pip install --no-cache-dir -r requirements-test.txt
COPY tests tests
RUN python -m pytest -q tests

FROM base AS run
RUN useradd --system --uid 10001 qsim
USER qsim
EXPOSE 8000
HEALTHCHECK --interval=15s --timeout=5s --start-period=10s --retries=3 \
  CMD python -c "import urllib.request,sys; sys.exit(0 if urllib.request.urlopen('http://127.0.0.1:8000/health').status==200 else 1)"
CMD ["uvicorn", "qsim.app:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "2", "--no-server-header"]
