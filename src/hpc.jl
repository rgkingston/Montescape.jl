# Generate a SLURM submission script for an HPC run.
#
# Key SLURM-layout principle for this pipeline:
#   - One task, many CPUs (--ntasks=1 --cpus-per-task=N).
#   - Parallelism (mc_workers, cs_threads) is set in config.yml by the
#     user and read at runtime. This script does NOT override those values.
#   - --cpus-per-task should be >= mc_workers * cs_threads. A warning is
#     emitted if the config implies more cores than requested here.
#   - This script only controls cluster resource booking:
#       cpus, mem_gb, walltime, partition, account.
"""
    write_slurm_script(config_path; cpus=64, mem_gb=128, walltime="24:00:00",
                       partition=nothing, account=nothing,
                       julia_exe="julia",
                       project_path=nothing, out_path=nothing) -> String

Write a .sbatch file next to the config and return its path.
Parallelism (mc_workers, cs_threads) is read directly from config.yml -
do not set it here. To change worker counts, edit the config in the
Montescape UI and regenerate this script.
The only parameters here are cluster resource requests:
  cpus      - --cpus-per-task (should be >= mc_workers * cs_threads)
  mem_gb    - memory in GB
  walltime  - wall clock limit (HH:MM:SS)
  partition - optional SLURM partition
  account   - optional SLURM account
"""
function write_slurm_script(config_path::AbstractString;
                             cpus::Integer = 64,
                             mem_gb::Integer = 128,
                             walltime::AbstractString = "24:00:00",
                             partition::Union{Nothing,AbstractString}=nothing,
                             account::Union{Nothing,AbstractString}=nothing,
                             julia_exe::AbstractString = "julia",
                             project_path::Union{Nothing,AbstractString}=nothing,
                             out_path::Union{Nothing,AbstractString}=nothing)

    cfg_path = abspath(config_path)
    cfg  = YAML.load_file(cfg_path; dicttype=Dict{String,Any})
    name = splitext(basename(cfg_path))[1]
    out  = out_path === nothing ? joinpath(dirname(cfg_path), "$name.sbatch") : out_path
    proj = project_path === nothing ? dirname(dirname(cfg_path)) : project_path

    # Read parallelism from the config - not overridden here
    cs_threads = Int(get(get(cfg, "circuitscape", Dict()), "julia_threads", 2))
    mc_workers = Int(get(get(cfg, "monte_carlo",  Dict()), "n_workers",     1))
    total = mc_workers * cs_threads

    if total > cpus
        @warn "SLURM script: mc_workers * cs_threads ($total) exceeds --cpus-per-task ($cpus)" mc_workers cs_threads cpus
    end

    open(out, "w") do io
        println(io, "#!/bin/bash")
        println(io, "# Montescape SLURM job - parallelism set in config.yml")
        println(io, "#SBATCH --job-name=$name")
        println(io, "#SBATCH --output=$(dirname(out))/$name.%j.out")
        println(io, "#SBATCH --error=$(dirname(out))/$name.%j.err")
        println(io, "#SBATCH --time=$walltime")
        println(io, "#SBATCH --ntasks=1")
        println(io, "#SBATCH --cpus-per-task=$cpus")
        println(io, "#SBATCH --mem=$(mem_gb)G")
        partition !== nothing && println(io, "#SBATCH --partition=$partition")
        account   !== nothing && println(io, "#SBATCH --account=$account")
        println(io)
        println(io, "set -euo pipefail")
        println(io)
        println(io, "# Parallelism is read from config.yml at runtime.")
        println(io, "# mc_workers=$mc_workers, cs_threads=$cs_threads, total=$total of $cpus cores")
        println(io, "# To change workers/threads edit config.yml, not this script.")
        println(io, "export JULIA_NUM_THREADS=$cs_threads")
        println(io, "export OPENBLAS_NUM_THREADS=$cs_threads")
        println(io)
        println(io, "echo \"==> Job started: \$(date)\"")
        println(io, "echo \"==> Host: \$(hostname)\"")
        println(io, "echo \"==> Config: $cfg_path\"")
        println(io, "echo \"==> Parallelism (from config): mc_workers=$mc_workers, cs_threads=$cs_threads, total=$total of $cpus cores\"")
        println(io)
        println(io, "cd $proj")
        println(io, "$julia_exe --project=. -t $cs_threads -e 'using MovementFramework; run_pipeline(\"$cfg_path\")'")
        println(io)
        println(io, "echo \"==> Job finished: \$(date)\"")
    end
    chmod(out, 0o755)
    return out
end