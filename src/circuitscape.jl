# Native Circuitscape.jl integration.
#

using Circuitscape
using LinearAlgebra

"""
    run_circuitscape(friction, ref, source_xy, ground_xy, run_dir; cfg, name) -> Dict

Applies circuitscape.julia_threads to the BLAS thread pool so the user's
UI / YAML setting actually controls the local Circuitscape parallelism.
Previously this field only affected the HPC sbatch and was silently
ignored on local runs.
"""
function run_circuitscape(friction::AbstractMatrix{<:Real},
                          ref::RasterData,
                          source_xy::Tuple{<:Real,<:Real},
                          ground_xy::Tuple{<:Real,<:Real},
                          run_dir::AbstractString;
                          cfg::AbstractDict,
                          name::AbstractString = "cs")

    cs_cfg = cfg["circuitscape"]
    # Honour circuitscape.julia_threads for local runs by pinning the BLAS
    # thread pool around this call. Reset on exit so this never leaks into
    # other parts of the pipeline.
    n_th = Int(get(cs_cfg, "julia_threads", 1))
    prev_blas = BLAS.get_num_threads()
    n_th > 0 && BLAS.set_num_threads(n_th)
    out_dir = joinpath(run_dir, "circuitscape", name)
    mkpath(out_dir)

    # Write resistance as ASCII grid (Circuitscape's preferred format).
    res_asc = joinpath(out_dir, "resistance.asc")
    _write_ascii_grid(res_asc, friction, ref)

    # Source and ground rasters as single-pixel masks.
    src_arr = zeros(Float32, size(friction))
    gnd_arr = zeros(Float32, size(friction))
    sr, sc = world_to_pixel(ref, source_xy...)
    gr, gc = world_to_pixel(ref, ground_xy...)
    src_arr[sr, sc] = 1.0f0
    gnd_arr[gr, gc] = 1.0f0
    src_asc = joinpath(out_dir, "source.asc")
    gnd_asc = joinpath(out_dir, "ground.asc")
    _write_ascii_grid(src_asc, src_arr, ref)
    _write_ascii_grid(gnd_asc, gnd_arr, ref)

    ini_path = joinpath(out_dir, "circuitscape.ini")
    _write_cs_ini(ini_path, cs_cfg, out_dir, res_asc, src_asc, gnd_asc, name)

    @info "  [circuitscape] computing $(name) with BLAS threads=$n_th..."
    try
        Circuitscape.compute(ini_path)
    finally
        BLAS.set_num_threads(prev_blas)
    end

    cur_map = joinpath(out_dir, "$(name)_curmap.asc")
    return Dict(
        "ini"     => ini_path,
        "out_dir" => out_dir,
        "curmap"  => isfile(cur_map) ? cur_map : "",
    )
end

function _write_ascii_grid(path::AbstractString, data::AbstractMatrix, ref::RasterData)
    rows, cols = size(data)
    gt = ref.geotransform
    xll = gt[1]
    yll = gt[4] + rows * gt[6]   # gt[6] is negative
    cellsize = abs(gt[2])
    open(path, "w") do io
        println(io, "ncols         $cols")
        println(io, "nrows         $rows")
        println(io, "xllcorner     $xll")
        println(io, "yllcorner     $yll")
        println(io, "cellsize      $cellsize")
        println(io, "NODATA_value  $(ref.nodata)")
        for r in 1:rows
            println(io, join((data[r, c] for c in 1:cols), " "))
        end
    end
end

function _write_cs_ini(path, cs_cfg, out_dir, res_path, src_path, gnd_path, name)
    scenario = get(cs_cfg, "scenario", "advanced")
    open(path, "w") do io
        println(io, "[Options for advanced mode]")
        println(io, "ground_file_is_resistances = False")
        println(io, "source_file = $src_path")
        println(io, "remove_src_or_gnd = keepall")
        println(io, "ground_file = $gnd_path")
        println(io, "use_unit_currents = $(get(cs_cfg, "use_unit_currents", true))")
        println(io, "use_direct_grounds = $(get(cs_cfg, "use_direct_grounds", true))")
        println(io)
        println(io, "[Calculation options]")
        println(io, "low_memory_mode = False")
        println(io, "parallelize = False")
        println(io, "solver = cg+amg")
        println(io, "print_timings = True")
        println(io)
        println(io, "[Short circuit regions (aka polygons)]")
        println(io, "use_polygons = False")
        println(io)
        println(io, "[Options for one-to-all and all-to-one modes]")
        println(io, "use_variable_source_strengths = False")
        println(io)
        println(io, "[Output options]")
        println(io, "set_null_currents_to_nodata = True")
        println(io, "set_null_voltages_to_nodata = True")
        println(io, "compress_grids = False")
        println(io, "write_cur_maps = $(get(cs_cfg, "write_cur_maps", true))")
        println(io, "write_volt_maps = $(get(cs_cfg, "write_volt_maps", false))")
        println(io, "output_file = $(joinpath(out_dir, "$(name).out"))")
        println(io, "write_cum_cur_map_only = False")
        println(io, "log_transform_maps = $(get(cs_cfg, "log_transform_maps", false))")
        println(io)
        println(io, "[Habitat raster or graph]")
        println(io, "habitat_file = $res_path")
        println(io, "habitat_map_is_resistances = True")
        println(io)
        println(io, "[Connection scheme for raster habitat data]")
        println(io, "connect_four_neighbors_only = $(get(cs_cfg, "connect_four_neighbors_only", false))")
        println(io, "connect_using_avg_resistances = $(get(cs_cfg, "use_average_resistances", true))")
        println(io)
        println(io, "[Version]")
        println(io, "version = 5.0.0")
        println(io)
        println(io, "[Mask file]")
        println(io, "use_mask = False")
        println(io)
        println(io, "[Circuitscape mode]")
        println(io, "data_type = raster")
        println(io, "scenario = $scenario")
    end
end
