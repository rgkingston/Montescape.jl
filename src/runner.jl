# Spawn the pipeline in a NEW detached terminal window.
#
# The GUI calls this when the user clicks "Run locally". The point is
# that the pipeline runs in its own console with live log output the
# user can watch, while the GUI process stays responsive and the user
# can configure / launch the next run.
#
# Windows quoting note:
#   `cmd /c start "title" cmd /k <command>` parses arguments in layers,
#   each of which strips one round of quotes. Embedding
#   `julia -e "using ..."` directly in that chain corrupts the inner
#   quotes and Julia sees `using` as a bare token, yielding
#   `ParseError: unterminated string literal`. We avoid the whole mess
#   by writing a tiny .bat file with the full command and spawning the
#   terminal on that file. Bash systems don't need this trick.

"""
    spawn_pipeline_in_terminal(config_path; project_dir=nothing,
                                julia_exe=joinpath(Sys.BINDIR, "julia"))

Launches a new terminal window that runs:
    julia --project=<dir> -e 'using MovementFramework; run_pipeline("<cfg>")'

Uses a .bat shim on Windows to avoid cmd.exe quoting hell. Uses an
available terminal emulator on Linux and `osascript` on macOS. Returns
nothing. Failures to find a terminal emulator on Linux fall back to
running detached without a visible window.
"""
function spawn_pipeline_in_terminal(config_path::AbstractString;
                                     project_dir::Union{Nothing,AbstractString}=nothing,
                                     julia_exe::AbstractString=joinpath(Sys.BINDIR, Sys.iswindows() ? "julia.exe" : "julia"))
    cfg = abspath(config_path)
    proj = project_dir === nothing ?
            dirname(dirname(@__DIR__)) :  # repo root, two up from src/
            abspath(project_dir)

    if Sys.iswindows()
        # Write the command to a .bat file in a temp dir, then spawn cmd /k on it.
        # cfg may contain backslashes; double them for the Julia string literal.
        cfg_for_julia = replace(cfg, "\\" => "\\\\")
        bat_path = joinpath(tempdir(), "mf_run_$(Base.Libc.getpid())_$(rand(UInt32)).bat")
        open(bat_path, "w") do io
            println(io, "@echo off")
            println(io, "title MovementFramework run - $(basename(cfg))")
            println(io, "echo Running pipeline on:")
            println(io, "echo   $(cfg)")
            println(io, "echo Using project:")
            println(io, "echo   $(proj)")
            println(io, "echo.")
            # Quote both paths defensively; parentheses in $proj (e.g. Downloads\foo(1)\)
            # would otherwise be interpreted by cmd.exe.
            println(io, "\"$(julia_exe)\" --project=\"$(proj)\" -e \"using MovementFramework; run_pipeline(\\\"$(cfg_for_julia)\\\")\"")
            println(io, "echo.")
            println(io, "echo Pipeline finished. Press any key to close this window.")
            println(io, "pause >nul")
        end
        @info "Spawning pipeline via bat shim" bat_path proj cfg
        run(Cmd(["cmd", "/c", "start", "", "cmd", "/k", bat_path]); wait=false)
        return nothing
    end

    # Unix-y systems: build the command string directly.
    julia_cmd = "$(julia_exe) --project=$(proj) -e 'using MovementFramework; run_pipeline(\"$(cfg)\")'"
    @info "spawning pipeline in new terminal" julia_cmd

    if Sys.isapple()
        applescript = """
            tell application "Terminal"
                activate
                do script "$(julia_cmd)"
            end tell
        """
        run(Cmd(["osascript", "-e", applescript]); wait=false)
    else
        # Linux: try common terminal emulators in order.
        for term in ("x-terminal-emulator", "gnome-terminal", "konsole",
                      "xfce4-terminal", "xterm")
            exe = Sys.which(term)
            exe === nothing && continue
            try
                if term == "gnome-terminal"
                    run(Cmd([exe, "--", "bash", "-c", "$(julia_cmd); exec bash"]); wait=false)
                else
                    run(Cmd([exe, "-e", "bash -c '$(julia_cmd); exec bash'"]); wait=false)
                end
                return nothing
            catch err
                @warn "failed to launch terminal" term err
            end
        end
        @warn "no terminal emulator found; running detached without console"
        run(Cmd(["bash", "-c", julia_cmd]); wait=false)
    end
    return nothing
end
