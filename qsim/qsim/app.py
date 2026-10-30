"""quanterall.com's quantum demo API on the NumPy engine (replaces the BEAM service behind /quantum/).

    GET /quantum/{bell,teleport,grover,qml,shor15}?shots=N
    GET /quantum/qec/{bitflip,phaseflip,shor9,steane7,surface3}?shots=N&error_qubit=Q|none
    GET /health

Same paths, parameters and answers as before: shots 1..10000 (anything else: 1024); at most four
simulations at a time (a fifth gets 503 at once); 30 s per simulation.
"""
from __future__ import annotations

import asyncio
import os

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from . import demos

DEFAULT_SHOTS = 1024
MAX_SHOTS = 10_000
MAX_CONCURRENT = int(os.environ.get("QSIM_MAX_CONCURRENT", "4"))
TIMEOUT_S = float(os.environ.get("QSIM_TIMEOUT_S", "30"))

CORS = {"access-control-allow-origin": "*",
        "access-control-allow-methods": "GET, OPTIONS",
        "access-control-allow-headers": "content-type"}

app = FastAPI(title="Quanterall quantum demos", docs_url=None, redoc_url=None, openapi_url=None)
_slots = asyncio.Semaphore(MAX_CONCURRENT)
_busy = 0


def _shots(raw) -> int:
    try:
        n = int(raw)
    except (TypeError, ValueError):
        return DEFAULT_SHOTS
    return n if 0 < n <= MAX_SHOTS else DEFAULT_SHOTS


def _error(status: int, message: str) -> JSONResponse:
    return JSONResponse({"status": "error", "message": message}, status_code=status, headers=CORS)


async def _simulate(example: str, fn, *args) -> JSONResponse:
    global _busy
    if _busy >= MAX_CONCURRENT:
        return _error(503, f"Server busy — max {MAX_CONCURRENT} concurrent simulations. Retry shortly.")
    _busy += 1
    try:
        async with _slots:
            result = await asyncio.wait_for(asyncio.to_thread(fn, *args), TIMEOUT_S)
        return JSONResponse(result | {"example": example}, headers=CORS)
    except asyncio.TimeoutError:
        return _error(503, f"Simulation timed out after {TIMEOUT_S:.0f}s.")
    except Exception:                                              # the details stay in the log
        import logging
        logging.getLogger("qsim").exception("simulation %s failed", example)
        return _error(500, "Internal simulation error.")
    finally:
        _busy -= 1


@app.get("/health")
async def health():
    return JSONResponse({"status": "ok", "engine": demos.ENGINE, "version": demos.VERSION, "busy": _busy},
                        headers=CORS)


@app.get("/quantum/qec/{name}")
async def qec(name: str, request: Request):
    fn = demos.QEC.get(name)
    if fn is None:
        return _error(404, "unknown demo")
    q = request.query_params
    return await _simulate(name, fn, _shots(q.get("shots")), q.get("error_qubit", "none"))


@app.get("/quantum/{name}")
async def ai(name: str, request: Request):
    fn = demos.AI.get(name)
    if fn is None:
        return _error(404, "unknown demo")
    return await _simulate(name, fn, _shots(request.query_params.get("shots")))


@app.options("/quantum/{path:path}")
async def preflight(path: str):
    return JSONResponse({}, headers=CORS)
