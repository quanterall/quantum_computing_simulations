-module(quanterall_quantum_router).
-behaviour(nova_router).
-export([routes/1]).

routes(_Environment) ->
    [
        #{
            prefix => "",
            security => false,
            routes => [
                {"/", {static_controller, index}, #{methods => [get]}},
                {"/health", {static_controller, health}, #{methods => [get]}}
            ]
        },
        #{
            prefix => "/quantum",
            security => false,
            routes => [
                {"/bell",         {quantum_controller, bell},      #{methods => [get]}},
                {"/teleport",     {quantum_controller, teleport},  #{methods => [get]}},
                {"/entangle",     {quantum_controller, entangle},  #{methods => [get]}},
                {"/migrate",      {quantum_controller, migrate},   #{methods => [get]}},
                {"/grover",       {quantum_controller, grover},    #{methods => [get]}},
                {"/qml",          {quantum_controller, qml},       #{methods => [get]}},
                {"/shor15",       {quantum_controller, shor15},    #{methods => [get]}},
                {"/optics/hom",           {quantum_controller, optics_hom},    #{methods => [get]}},
                {"/optics/ns",            {quantum_controller, optics_ns},     #{methods => [get]}},
                {"/optics/fusion",        {quantum_controller, optics_fusion}, #{methods => [get]}},
                {"/optics/cz",            {quantum_controller, optics_cz},     #{methods => [get]}},
                {"/optics/gkp",           {quantum_controller, optics_gkp},    #{methods => [get]}},
                {"/optics/cat",           {quantum_controller, optics_cat},    #{methods => [get]}},
                {"/optics/comb",          {quantum_controller, optics_comb},   #{methods => [get]}},
                {"/shor",                 {quantum_controller, shor_general},  #{methods => [get]}},
                {"/grover/search",        {quantum_controller, grover_search}, #{methods => [get]}},
                {"/stabilizer/ghz",       {quantum_controller, stab_ghz},      #{methods => [get]}},
                {"/stabilizer/memory",    {quantum_controller, stab_memory},   #{methods => [get]}},
                {"/qec/bitflip",  {quantum_controller, bitflip},   #{methods => [get]}},
                {"/qec/phaseflip",{quantum_controller, phaseflip}, #{methods => [get]}},
                {"/qec/shor9",    {quantum_controller, shor9},     #{methods => [get]}},
                {"/qec/steane7",  {quantum_controller, steane7},   #{methods => [get]}},
                {"/qec/surface3", {quantum_controller, surface3},  #{methods => [get]}}
            ]
        }
    ].
