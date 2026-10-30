-module(quantum_controller).
-export([bell/1, teleport/1, grover/1, qml/1, shor15/1, entangle/1, migrate/1,
         grover_search/1, stab_ghz/1, stab_memory/1, shor_general/1, optics_hom/1, optics_ns/1, optics_cz/1, optics_fusion/1, optics_gkp/1, optics_cat/1, optics_comb/1,
         bitflip/1, phaseflip/1, shor9/1, steane7/1, surface3/1]).

-define(DEFAULT_SHOTS, 1024).
-define(TIMEOUT_MS,    30000).  % 30 seconds max per simulation

%% ── AI Demos ──────────────────────────────────────────────────────────────────
bell(Req)     -> call_demo(Req, bell,     fun 'Elixir.QxDemo':run_bell/1).
teleport(Req) -> call_demo(Req, teleport, fun 'Elixir.QxDemo':run_teleport/1).
grover(Req)   -> call_demo(Req, grover,   fun 'Elixir.QxDemo':run_grover/1).
qml(Req)      -> call_demo(Req, qml,      fun 'Elixir.QxDemo':run_qml/1).
shor15(Req)   -> call_demo(Req, shor15,   fun 'Elixir.QxDemo':run_shor15/1).

%% ── Grover on n qubits; the stabilizer simulator (lib/qx_stab.ex). Every parameter is bounded. ──
grover_search(Req) ->
    P = qs(Req),
    N = int_param(P, <<"n">>, 4, 2, 10),
    Marked = int_param(P, <<"marked">>, (1 bsl N) - 1, 0, (1 bsl N) - 1),
    Iter = case maps:get(<<"iterations">>, P, undefined) of
               undefined -> auto;
               _ -> int_param(P, <<"iterations">>, 1, 0, 80)
           end,
    call_demo(Req, grover_search, fun(S) -> 'Elixir.QxDemo':run_grover_search(N, Marked, Iter, S) end).

%% Photonic schemes from the literature (lib/qx_fock.ex): HOM, the KLM NS gate and CZ
optics_hom(Req) ->
    P = qs(Req),
    GT = float_param(P, <<"gt">>, 0.0, -10.0, 10.0),
    call_demo(Req, optics_hom, fun(_S) -> 'Elixir.QxFock':run_hom(GT) end).

optics_ns(Req) ->
    call_demo(Req, optics_ns, fun(_S) -> 'Elixir.QxFock':run_ns() end).

optics_fusion(Req) ->
    call_demo(Req, optics_fusion, fun(S) -> 'Elixir.QxFock':run_fusion(min(S, 20000)) end).

optics_cz(Req) ->
    P = qs(Req),
    In = case maps:get(<<"input">>, P, <<"11">>) of I when I =:= <<"00">>; I =:= <<"01">>; I =:= <<"10">>; I =:= <<"11">> -> I; _ -> <<"11">> end,
    Off = float_param(P, <<"theta2_offset">>, 0.0, -5.0, 5.0),
    Phase = maps:get(<<"phase">>, P, <<"1">>) =/= <<"0">>,
    Opts = [{angles, {22.5, 65.5302 + Off, -22.5}}, {phase_pi, Phase}],
    call_demo(Req, optics_cz, fun(S) -> 'Elixir.QxFock':run_cz(In, min(S, 20000), Opts) end).

%% Phase 7, one mode (lib/qx_cv.ex): GKP readout under shift noise, cat parity under loss, the heralded comb
optics_gkp(Req) ->
    P = qs(Req),
    Db = float_param(P, <<"db">>, 10.0, 3.0, 12.0),
    D2 = float_param(P, <<"delta2">>, 0.01, 0.005, 0.2),
    Sigma = float_param(P, <<"sigma">>, 0.3, 0.0, 1.0),
    Basis = case maps:get(<<"basis">>, P, <<"z">>) of <<"x">> -> x; _ -> z end,
    call_demo(Req, optics_gkp, fun(S) -> 'Elixir.QxCV':run_gkp(Db, D2, Sigma, Basis, min(S, 20000)) end).

optics_cat(Req) ->
    P = qs(Req),
    Beta = float_param(P, <<"beta">>, 2.0, 0.3, 4.0),
    Parity = case maps:get(<<"parity">>, P, <<"even">>) of <<"odd">> -> 1; _ -> 0 end,
    Loss = float_param(P, <<"loss">>, 0.05, 0.0, 0.5),
    call_demo(Req, optics_cat, fun(S) -> 'Elixir.QxCV':run_cat(Beta, Parity, 1.0 - Loss, min(S, 20000)) end).

optics_comb(Req) ->
    P = qs(Req),
    R = float_param(P, <<"r">>, 1.5, 0.5, 2.0),
    Beta = float_param(P, <<"beta">>, 3.0, 1.0, 5.0),
    Theta = float_param(P, <<"theta">>, math:sqrt(math:pi()), 0.5, 3.0),
    Phi = float_param(P, <<"phi">>, 0.0, -math:pi(), math:pi()),
    Y = case maps:is_key(<<"y">>, P) of true -> float_param(P, <<"y">>, 0.0, -20.0, 20.0); false -> nil end,
    call_demo(Req, optics_comb, fun(_S) -> 'Elixir.QxCV':run_comb(R, Beta, Theta, Phi, Y) end).

shor_general(Req) ->
    P = qs(Req),
    N = case int_param(P, <<"N">>, 15, 15, 21) of 21 -> 21; _ -> 15 end,
    A = int_param(P, <<"a">>, 7, 2, N - 1),                % a sharing a factor with N is answered at once
    call_demo(Req, shor, fun(S) -> 'Elixir.QxDemo':run_shor(N, A, S) end).

gcd(A, 0) -> A;
gcd(A, B) -> gcd(B, A rem B).

stab_ghz(Req) ->
    P = qs(Req),
    N = int_param(P, <<"n">>, 1000, 2, 2000),
    call_demo(Req, stabilizer_ghz, fun(S) -> 'Elixir.QxStab':run_ghz(N, S) end).

stab_memory(Req) ->
    P = qs(Req),
    Code = case maps:get(<<"code">>, P, <<"surface">>) of <<"repetition">> -> repetition; _ -> surface end,
    D = int_param(P, <<"d">>, 3, 3, 7),
    Rounds = int_param(P, <<"rounds">>, D, 1, 3 * D),
    Prob = float_param(P, <<"p">>, 0.001, 0.0, 0.05),
    call_demo(Req, stabilizer_memory, fun(S) -> 'Elixir.QxStab':run_memory(Code, D, Rounds, Prob, S) end).

qs(Req) -> maps:from_list(cowboy_req:parse_qs(Req)).

int_param(P, Key, Default, Min, Max) ->
    case safe_to_integer(maps:get(Key, P, <<>>)) of
        {ok, V} when V >= Min, V =< Max -> V;
        _ -> Default
    end.

float_param(P, Key, Default, Min, Max) ->
    Bin = maps:get(Key, P, <<>>),
    V = try binary_to_float(Bin) catch _:_ -> try float(binary_to_integer(Bin)) catch _:_ -> Default end end,
    if V >= Min, V =< Max -> V; true -> Default end.

%% ── Process models (lib/qx_processes.ex): "shots" is the number of rounds / runs ──
entangle(Req) -> call_demo(Req, entangle, fun(N) -> 'Elixir.QxProcesses':run_entangle(min(N, 2000)) end).
migrate(Req)  -> call_demo(Req, migrate,  fun(N) -> 'Elixir.QxProcesses':run_teleport(min(N, 300)) end).

%% ── QEC Demos ─────────────────────────────────────────────────────────────────
bitflip(Req)   -> call_qec(Req, bitflip,   fun 'Elixir.QxDemo':run_bitflip_repetition/2).
phaseflip(Req) -> call_qec(Req, phaseflip, fun 'Elixir.QxDemo':run_phaseflip_repetition/2).
shor9(Req)     -> call_qec(Req, shor9,     fun 'Elixir.QxDemo':run_shor9/2).
steane7(Req)   -> call_qec(Req, steane7,   fun 'Elixir.QxDemo':run_steane7/2).
surface3(Req)  -> call_qec(Req, surface3,  fun 'Elixir.QxDemo':run_surface3/2).

%% ── Internal: AI demo dispatcher ──────────────────────────────────────────────
call_demo(Req, Example, Fun) ->
    Params = maps:from_list(cowboy_req:parse_qs(Req)),
    Shots  = parse_shots(maps:get(<<"shots">>, Params, undefined)),
    Seed   = parse_seed(maps:get(<<"seed">>, Params, undefined)),
    Headers = cors_headers(),
    %% Fix 2: acquire a concurrency slot before running
    case 'Elixir.QuantumLimiter':acquire() of
        ok ->
            Result = run_with_timeout(Seed, fun() -> Fun(Shots) end),
            'Elixir.QuantumLimiter':release(),
            case Result of
                {ok, R}  -> {json, 200, Headers, R#{example => Example, seed => Seed}};
                timeout  -> {json, 503, Headers, busy_error(<<"Simulation timed out after 30s.">>)};
                {error,_}-> {json, 500, Headers, busy_error(<<"Internal simulation error.">>)}
            end;
        {error, busy} ->
            {json, 503, Headers, busy_error(<<"Server busy — max 4 concurrent simulations. Retry shortly.">>)};
        {error, _} ->
            {json, 500, Headers, busy_error(<<"Limiter unavailable.">>)}
    end.

%% ── Internal: QEC demo dispatcher ─────────────────────────────────────────────
call_qec(Req, Example, Fun) ->
    Params     = maps:from_list(cowboy_req:parse_qs(Req)),
    Shots      = parse_shots(maps:get(<<"shots">>, Params, undefined)),
    ErrorQubit = maps:get(<<"error_qubit">>, Params, <<"none">>),
    Seed       = parse_seed(maps:get(<<"seed">>, Params, undefined)),
    Headers    = cors_headers(),
    case 'Elixir.QuantumLimiter':acquire() of
        ok ->
            Result = run_with_timeout(Seed, fun() -> Fun(Shots, ErrorQubit) end),
            'Elixir.QuantumLimiter':release(),
            case Result of
                {ok, R}  -> {json, 200, Headers, R#{example => Example, seed => Seed}};
                timeout  -> {json, 503, Headers, busy_error(<<"QEC simulation timed out after 30s.">>)};
                {error,_}-> {json, 500, Headers, busy_error(<<"Internal simulation error.">>)}
            end;
        {error, busy} ->
            {json, 503, Headers, busy_error(<<"Server busy — max 4 concurrent simulations. Retry shortly.">>)};
        {error, _} ->
            {json, 500, Headers, busy_error(<<"Limiter unavailable.">>)}
    end.

%% ── Fix 2: run simulation in a separate process with hard timeout ──────────────
%% The simulation's random numbers come from Seed, so the same seed, shots and parameters give the same counts.
%% (The process models spawn many processes of their own and are not covered by it.)
run_with_timeout(Seed, Fun) ->
    Caller = self(),
    Ref    = make_ref(),
    Pid    = spawn(fun() ->
        try
            _ = rand:seed(exsss, Seed),
            Caller ! {Ref, ok, Fun()}
        catch
            _Class:Reason ->
                Caller ! {Ref, error, Reason}
        end
    end),
    receive
        {Ref, ok,    Result} -> {ok, Result};
        {Ref, error, Reason} -> {error, Reason}
    after ?TIMEOUT_MS ->
        exit(Pid, kill),
        %% Drain any late message to avoid mailbox pollution
        receive {Ref, _, _} -> ok after 0 -> ok end,
        timeout
    end.

%% ── Helpers ───────────────────────────────────────────────────────────────────
parse_shots(undefined) -> ?DEFAULT_SHOTS;
parse_shots(Bin) when is_binary(Bin) ->
    case safe_to_integer(Bin) of
        {ok, N} when N > 0, N =< 10000 -> N;
        _                               -> ?DEFAULT_SHOTS
    end.

%% a seed from the query (0 .. 2^32-1), or a fresh one; either way the response reports it
parse_seed(Bin) when is_binary(Bin) ->
    case safe_to_integer(Bin) of
        {ok, N} when N >= 0, N < 4294967296 -> N;
        _                                   -> parse_seed(undefined)
    end;
parse_seed(_) -> rand:uniform(4294967296) - 1.

safe_to_integer(Bin) ->
    try {ok, binary_to_integer(Bin)}
    catch _:_ -> error
    end.

cors_headers() ->
    #{
        <<"access-control-allow-origin">>  => <<"*">>,
        <<"access-control-allow-methods">> => <<"GET, OPTIONS">>,
        <<"access-control-allow-headers">> => <<"content-type">>
    }.

busy_error(Msg) ->
    #{status => <<"error">>, message => Msg}.
