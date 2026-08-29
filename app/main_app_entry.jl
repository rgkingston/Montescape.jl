# Entry-point shim referenced by PackageCompiler in build/build_app.jl.
# When the standalone binary launches, this starts the web GUI server
# and opens the user's default browser.

module MovementFrameworkApp

function julia_main()::Cint
    try
        include(joinpath(@__DIR__, "desktop_web.jl"))
        Main.main()
        return 0
    catch err
        @error "MovementFramework crashed" exception = (err, catch_backtrace())
        return 1
    end
end

end # module

using .MovementFrameworkApp: julia_main
