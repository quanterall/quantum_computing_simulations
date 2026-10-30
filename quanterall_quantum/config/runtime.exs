import Config

config :nova,
  bootstrap_application: :quanterall_quantum,
  router: :quanterall_quantum_router,
  environment: :prod,
  dev_mode: false,
  cowboy_configuration: %{
    port: String.to_integer(System.get_env("PORT") || "8080")
  }
