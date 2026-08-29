# Self-contained run package generator.
#
# Produces a folder the user can copy to an HPC cluster, with:
#   - the YAML config (paths rewritten relative to the package root)
#   - any custom_layer_files referenced from the config
#   - the package's source code (MovementFramework + Project.toml)
#   - the input rasters/vectors referenced in data.{dem_raster,
#     sites_vector, passes_vector}, copied into ./inputs/
#   - submit_slurm.sh, run_local.sh, watch_job.sh, README_HPC.txt

using YAML
using Dates

"""
    generate_run_package(config_path; output_dir=".", kwargs...) -> String

Returns the path to the generated package directory.
"""
function generate_run_package(config_path::AbstractString;
                               output_dir::AbstractString = ".",
                               package_name::Union{Nothing,AbstractString} = nothing,
                               cpus::Integer = 64,
                               mem_gb::Integer = 128,
                               walltime::AbstractString = "24:00:00",
                               partition::Union{Nothing,AbstractString} = nothing,
                               account::Union{Nothing,AbstractString} = nothing,
                               cs_threads::Integer = 2,
                               mc_workers::Integer = 32,
                               julia_module::AbstractString = "julia",
                               scratch_base::AbstractString = "/scratch/\$USER",
                               hpc_user::AbstractString = "yourname",
                               hpc_host::AbstractString = "hpc.cluster.example",
                               include_inputs::Bool = true,
                               julia_threads::Union{Nothing,Integer}=nothing)

    if julia_threads !== nothing
        cs_threads = julia_threads
    end

    cfg_path = abspath(config_path)
    cfg = YAML.load_file(cfg_path; dicttype=Dict{String,Any})

    name = package_name === nothing ?
        string(cfg["project"]["name"], "_", Dates.format(now(), "yyyymmdd_HHMMSS")) :
        package_name
    pkg_dir = abspath(joinpath(output_dir, name))
    mkpath(pkg_dir)
    mkpath(joinpath(pkg_dir, "inputs"))
    mkpath(joinpath(pkg_dir, "custom_layers"))

    repo_root = _find_repo_root(cfg_path)

    if include_inputs
        for key in ("dem_raster", "sites_vector", "passes_vector")
            src = get(cfg["data"], key, nothing)
            src === nothing && continue
            src_abs = isabspath(src) ? src : abspath(joinpath(repo_root, src))
            isfile(src_abs) || (@warn "Skipping missing input" key src_abs; continue)
            dst_rel = joinpath("inputs", basename(src_abs))
            cp(src_abs, joinpath(pkg_dir, dst_rel); force=true)
            cfg["data"][key] = dst_rel
            if endswith(lowercase(src_abs), ".shp")
                base = first(splitext(src_abs))
                for ext in (".shx", ".dbf", ".prj", ".cpg")
                    side = base * ext
                    isfile(side) && cp(side, joinpath(pkg_dir, "inputs", basename(side)); force=true)
                end
            end
        end
    end

    new_custom_files = String[]
    for clf in get(cfg["friction"], "custom_layer_files", String[])
        clf_abs = isabspath(clf) ? clf : abspath(joinpath(repo_root, clf))
        if !isfile(clf_abs)
            @warn "Skipping missing custom layer file" clf clf_abs
            continue
        end
        dst_rel = joinpath("custom_layers", basename(clf_abs))
        cp(clf_abs, joinpath(pkg_dir, dst_rel); force=true)
        push!(new_custom_files, dst_rel)
    end
    cfg["friction"]["custom_layer_files"] = new_custom_files

    src_dir = joinpath(repo_root, "src")
    isdir(src_dir) && cp(src_dir, joinpath(pkg_dir, "src"); force=true)
    proj_toml = joinpath(repo_root, "Project.toml")
    isfile(proj_toml) && cp(proj_toml, joinpath(pkg_dir, "Project.toml"); force=true)

    # Override Monte Carlo workers + output_dir for HPC
    if haskey(cfg, "monte_carlo") && cfg["monte_carlo"] isa AbstractDict
        cfg["monte_carlo"]["n_workers"] = Int(mc_workers)
    end
    cfg["project"]["output_dir"] = joinpath(scratch_base, "outputs")

    cfg_out = joinpath(pkg_dir, "config.yml")
    YAML.write_file(cfg_out, cfg)

    _write_slurm_script(pkg_dir, name; cpus, mem_gb, walltime,
                        partition, account, cs_threads, mc_workers,
                        julia_module, scratch_base)
    _write_local_script(pkg_dir, name; cs_threads)
    _write_watch_job_script(pkg_dir, name; hpc_user, hpc_host)
    _write_hpc_readme(pkg_dir, name; hpc_user, hpc_host, scratch_base,
                      cpus, mem_gb, walltime, mc_workers, cs_threads)

    @info "Wrote run package: $pkg_dir"
    return pkg_dir
end

function _find_repo_root(cfg_path::AbstractString)
    d = dirname(abspath(cfg_path))
    for _ in 1:8
        isfile(joinpath(d, "Project.toml")) && return d
        parent = dirname(d)
        parent == d && break
        d = parent
    end
    return dirname(abspath(cfg_path))
end

function _write_slurm_script(dir, name; cpus, mem_gb, walltime,
                              partition, account, cs_threads, mc_workers,
                              julia_module, scratch_base)
    path = joinpath(dir, "submit_slurm.sh")
    total = Int(mc_workers) * Int(cs_threads)
    open(path, "w") do io
        println(io, "#!/bin/bash")
        println(io, "# Montescape SLURM job")
        println(io, "#SBATCH --job-name=$name")
        println(io, "#SBATCH --output=logs/$name.%j.out")
        println(io, "#SBATCH --error=logs/$name.%j.err")
        println(io, "#SBATCH --time=$walltime")
        println(io, "#SBATCH --ntasks=1")
        println(io, "#SBATCH --cpus-per-task=$cpus")
        println(io, "#SBATCH --mem=$(mem_gb)G")
        partition !== nothing && println(io, "#SBATCH --partition=$partition")
        account   !== nothing && println(io, "#SBATCH --account=$account")
        println(io)
        println(io, "# Edit this block if your cluster's scratch is NOT $scratch_base")
        println(io, "# (e.g. some clusters use /work/\$USER or /lustre/\$USER).")
        println(io, "SCRATCH_BASE=\"$scratch_base\"")
        println(io, "RUN_DIR=\"\${SCRATCH_BASE}/runs/$name\"")
        println(io, "OUTPUT_DIR=\"\${SCRATCH_BASE}/outputs\"")
        println(io)
        println(io, "set -euo pipefail")
        println(io, "mkdir -p \"\$RUN_DIR\" \"\$OUTPUT_DIR\" logs")
        println(io)
        println(io, "rsync -a --exclude logs --exclude '*.out' --exclude '*.err' \\")
        println(io, "      \"\$SLURM_SUBMIT_DIR/\" \"\$RUN_DIR/\"")
        println(io, "cd \"\$RUN_DIR\"")
        println(io)
        println(io, "module load $julia_module 2>/dev/null || true")
        println(io)
        println(io, "# Per-process thread count (BLAS / Circuitscape).")
        println(io, "# Monte Carlo workers are spawned from Julia.")
        println(io, "# Cores used: mc_workers * JULIA_NUM_THREADS = $mc_workers * $cs_threads = $total.")
        println(io, "export JULIA_NUM_THREADS=$cs_threads")
        println(io, "export OPENBLAS_NUM_THREADS=$cs_threads")
        println(io)
        println(io, "echo \"==> Job started: \$(date)\"")
        println(io, "echo \"==> Host: \$(hostname)\"")
        println(io, "echo \"==> Run dir: \$RUN_DIR\"")
        println(io, "echo \"==> Output dir: \$OUTPUT_DIR\"")
        println(io, "echo \"==> Parallelism: mc_workers=$mc_workers, cs_threads=$cs_threads, total=$total of $cpus cores\"")
        println(io)
        println(io, "julia --project=. -e 'using Pkg; Pkg.instantiate()'")
        println(io, "julia --project=. -t $cs_threads -e 'using MovementFramework; run_pipeline(\"config.yml\")'")
        println(io)
        println(io, "echo \"==> Job finished: \$(date)\"")
    end
    chmod(path, 0o755)
end

function _write_local_script(dir, name; cs_threads)
    path = joinpath(dir, "run_local.sh")
    open(path, "w") do io
        println(io, "#!/bin/bash")
        println(io, "# Sanity-check run on a local machine before sending to HPC.")
        println(io, "set -euo pipefail")
        println(io, "cd \"\$(dirname \"\$0\")\"")
        println(io, "export JULIA_NUM_THREADS=$cs_threads")
        println(io, "julia --project=. -e 'using Pkg; Pkg.instantiate()'")
        println(io, "julia --project=. -t $cs_threads -e 'using MovementFramework; run_pipeline(\"config.yml\")'")
    end
    chmod(path, 0o755)
end

# Helper that ssh's into the cluster and tails the latest job log so
# the user can see live output from a SLURM-submitted job in real time.
function _write_watch_job_script(dir, name; hpc_user, hpc_host)
    path = joinpath(dir, "watch_job.sh")
    open(path, "w") do io
        println(io, "#!/bin/bash")
        println(io, "# Watch a running (or finished) Montescape SLURM job in real time.")
        println(io, "#")
        println(io, "# After you sbatch on the cluster, run this from your laptop to")
        println(io, "# stream the job's stdout/stderr as it is written.")
        println(io, "#")
        println(io, "# Override the defaults by exporting these env vars:")
        println(io, "#   MS_HPC_USER, MS_HPC_HOST, MS_REMOTE_DIR")
        println(io)
        println(io, "USER_REMOTE=\"\${MS_HPC_USER:-$hpc_user}\"")
        println(io, "HOST_REMOTE=\"\${MS_HPC_HOST:-$hpc_host}\"")
        println(io, "REMOTE_DIR=\"\${MS_REMOTE_DIR:-~/runs/$name}\"")
        println(io)
        println(io, "echo \"Watching \$USER_REMOTE@\$HOST_REMOTE:\$REMOTE_DIR/logs/ (Ctrl+C to stop)\"")
        println(io, "ssh -t \"\$USER_REMOTE@\$HOST_REMOTE\" \"cd \$REMOTE_DIR && tail -F logs/*.out logs/*.err 2>/dev/null\"")
    end
    chmod(path, 0o755)
end

function _write_hpc_readme(dir, name; hpc_user, hpc_host, scratch_base,
                            cpus, mem_gb, walltime, mc_workers, cs_threads)
    open(joinpath(dir, "README_HPC.txt"), "w") do io
        println(io, "=========================================================")
        println(io, " Montescape run package: $name")
        println(io, "=========================================================")
        println(io)
        println(io, "This folder contains everything needed to run the")
        println(io, "analysis on an HPC cluster. To submit it:")
        println(io)
        println(io, "  # 1. From your local machine, copy the folder to the cluster:")
        println(io, "  scp -r \"\$(pwd)\" $hpc_user@$hpc_host:~/runs/")
        println(io)
        println(io, "  # 2. Log in and submit:")
        println(io, "  ssh $hpc_user@$hpc_host")
        println(io, "  cd ~/runs/$name")
        println(io, "  sbatch submit_slurm.sh")
        println(io)
        println(io, "  # 3. From your laptop, watch the job live:")
        println(io, "  ./watch_job.sh")
        println(io)
        println(io, "  # Or, on the cluster:")
        println(io, "  squeue -u \$USER")
        println(io, "  tail -F logs/$name.*.out")
        println(io)
        println(io, "Outputs will be written to: $scratch_base/outputs/<run_id>/")
        println(io)
        println(io, "---------------------------------------------------------")
        println(io, " Parallelism layout")
        println(io, "---------------------------------------------------------")
        println(io, "  mc_workers : $mc_workers   (Monte Carlo iterations in flight)")
        println(io, "  cs_threads : $cs_threads   (BLAS / Circuitscape threads per iter)")
        println(io, "  total cores: $(Int(mc_workers)*Int(cs_threads)) of --cpus-per-task=$cpus")
        println(io)
        println(io, "If your cluster reports low CPU usage in seff, lower")
        println(io, "cs_threads or raise mc_workers. They are independent dials.")
        println(io)
        println(io, "  CPUs per task : $cpus")
        println(io, "  Memory        : $(mem_gb) GB")
        println(io, "  Walltime      : $walltime")
        println(io)
        println(io, "  sbatch submit_slurm.sh        - submit the job")
        println(io, "  squeue -u \$USER               - show your jobs")
        println(io, "  scancel <jobid>               - cancel a job")
        println(io, "  seff <jobid>                  - efficiency after completion")
    end
end
