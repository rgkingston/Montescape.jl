# Minimal CLI entry point: julia --project=. src/cli.jl <config.yml>

function main_cli(args::AbstractVector{<:AbstractString}=ARGS)
    if isempty(args)
        println("Usage: julia --project=. -e 'using MovementFramework; run_pipeline(\"path/to/config.yml\")'")
        println("   or: julia --project=. src/cli.jl path/to/config.yml")
        return 1
    end
    run_pipeline(args[1])
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_cli(ARGS))
end
