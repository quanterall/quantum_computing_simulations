-module(static_controller).
-export([index/1, health/1]).

index(Req) ->
    FilePath = code:priv_dir(quanterall_quantum),
    FullPath = filename:join([FilePath, "static", "index.html"]),
    case file:read_file(FullPath) of
        {ok, Bin} ->
            {ok, Req2} = cowboy_req:reply(200,
                #{<<"content-type">> => <<"text/html; charset=utf-8">>},
                Bin, Req),
            {ok, Req2};
        _ ->
            {json, 404, #{}, #{error => <<"not_found">>}}
    end.

health(_Req) ->
    {json, 200, #{}, #{status => <<"ok">>, node => node()}}.
