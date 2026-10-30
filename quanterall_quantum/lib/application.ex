defmodule QuanterallQuantum.Application do
  use Application

  @impl true
  def start(_type, _args) do
    children = [
      # Fix 2: concurrency limiter — must start before the HTTP server
      {QuantumLimiter, []},
      # node B of the teleportation model: a second BEAM node (OTP :peer), started on first use
      {QxProcesses.NodeB, []}
    ]
    opts = [strategy: :one_for_one, name: QuanterallQuantum.Supervisor]
    Supervisor.start_link(children, opts)
  end
end
