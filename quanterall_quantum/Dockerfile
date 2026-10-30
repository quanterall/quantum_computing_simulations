# Elixir 1.18 on Erlang/OTP 27
FROM elixir:1.18

ENV MIX_ENV=prod
ENV PORT=8080

WORKDIR /app

RUN apt-get update && \
    apt-get install -y build-essential git curl && \
    rm -rf /var/lib/apt/lists/*

RUN mix local.hex --force && \
    mix local.rebar --force

# Resolve and compile deps
COPY mix.exs ./
RUN mix deps.get --only prod && \
    mix deps.compile

# Copy application source
COPY config config
COPY lib lib
COPY src src
COPY priv priv

# Compile the application (but NOT a release — nova needs the full OTP
# environment at runtime to dynamically compile nova.Router via compile:forms/2)
RUN mix compile

EXPOSE 8080

HEALTHCHECK --interval=15s --timeout=5s --start-period=15s --retries=3 \
  CMD curl -f http://localhost:8080/health || exit 1

CMD ["mix", "run", "--no-halt"]
