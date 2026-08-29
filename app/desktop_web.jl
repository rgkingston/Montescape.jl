# dev launcher for the MovementFramework web GUI.
#
# In a checked-out source tree, run:
#   julia --project=. app/desktop_web.jl
#

using MovementFramework

if abspath(PROGRAM_FILE) == @__FILE__
    MovementFramework.WebGUI.main()
end
