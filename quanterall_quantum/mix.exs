defmodule QuanterallQuantum.MixProject do
  use Mix.Project

  def project do
    [
      app: :quanterall_quantum,
      version: "3.0.0",
      elixir: "~> 1.18",
      deps: deps()
    ]
  end

  def application do
    [
      mod: {QuanterallQuantum.Application, []},
      extra_applications: [:logger]
    ]
  end

  defp deps do
    [
      {:nova, "~> 0.13.0"}
    ]
  end
end
