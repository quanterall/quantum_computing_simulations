%% Bootstrap stub for nova_router:compile/3.
%%
%% Nova 0.13.x calls 'nova.Router':module_info(compile) at startup
%% before generating the real routing module, assuming this module
%% already exists (e.g. from a prior hot-reload). When starting fresh
%% (as happens in Docker), the module doesn't exist and the call
%% crashes with undef.
%%
%% Compiling this stub gives Nova a live module to call module_info/1
%% on. nova_router then proceeds to compile the real routes and
%% replaces this stub in memory via code:load_binary/3.
%%
%% Do NOT export any functions — all Erlang modules automatically
%% export module_info/0 and module_info/1, which is all Nova needs here.
-module('nova.Router').
