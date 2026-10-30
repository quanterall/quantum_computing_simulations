%% The router Nova uses for this application. Nova 0.13 detects an Elixir (mix) project and then looks
%% for the routes in '<bootstrap_application>.Router' -- this module -- not in <app>_router. It used to
%% be an empty stub, so no route was registered and every path, /health included, answered 404.
%% The routes stay in quanterall_quantum_router; this module hands them over.
-module('quanterall_quantum.Router').
-export([routes/1]).

routes(Environment) ->
    quanterall_quantum_router:routes(Environment).
