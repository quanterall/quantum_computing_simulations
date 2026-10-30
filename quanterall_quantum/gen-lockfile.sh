#!/usr/bin/env bash
# Run once before building the image to generate/refresh mix.lock
# Uses the same Elixir image as the Dockerfile to guarantee compatibility.
set -euo pipefail
cd "$(dirname "$0")"

echo "Generating mix.lock using elixir:1.18..."
docker run --rm \
  -v "$(pwd):/app" \
  -w /app \
  elixir:1.18 \
  sh -c "mix local.hex --force --quiet && mix local.rebar --force --quiet && MIX_ENV=prod mix deps.get --only prod"

echo "Done — mix.lock written."
