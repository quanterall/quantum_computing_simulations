defmodule QuantumLimiter do
  @moduledoc """
  Fix 2: Concurrency limiter for quantum simulations.

  Caps simultaneous simulations at @max_concurrent (default 4).
  Uses a GenServer as the single source of truth — no race conditions.
  Erlang controller calls acquire/0 and release/0 via BEAM cross-language calls.
  """
  use GenServer

  @max_concurrent 4

  # ── Public API ──────────────────────────────────────────────────────────────

  def start_link(_opts) do
    GenServer.start_link(__MODULE__, 0, name: :quantum_limiter)
  end

  @doc "Try to acquire a simulation slot. Returns :ok or {:error, :busy}."
  def acquire do
    try do
      GenServer.call(:quantum_limiter, :acquire, 5_000)
    catch
      :exit, _ -> {:error, :unavailable}
    end
  end

  @doc "Release a previously acquired simulation slot."
  def release do
    GenServer.cast(:quantum_limiter, :release)
  end

  @doc "Return current concurrency count (for /health endpoint)."
  def count do
    try do
      GenServer.call(:quantum_limiter, :count, 2_000)
    catch
      :exit, _ -> -1
    end
  end

  # ── GenServer callbacks ──────────────────────────────────────────────────────

  @impl true
  def init(_), do: {:ok, 0}

  @impl true
  def handle_call(:acquire, _from, count) when count < @max_concurrent do
    {:reply, :ok, count + 1}
  end

  def handle_call(:acquire, _from, count) do
    {:reply, {:error, :busy}, count}
  end

  def handle_call(:count, _from, count) do
    {:reply, count, count}
  end

  @impl true
  def handle_cast(:release, count) do
    {:noreply, max(0, count - 1)}
  end
end
