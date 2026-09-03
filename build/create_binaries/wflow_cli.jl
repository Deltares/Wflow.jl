using wflow_cli

function (@main)(_)::Cint
    return wflow_cli.julia_main()
end
