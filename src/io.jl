# Raster + vector IO via ArchGDAL.
#
# IMPORTANT
# indexing convention:
#   ArchGDAL.read(band) returns an array indexed as arr[x, y] (i.e. arr[col, row])
#   following GDAL's column-major X/Y convention. The rest of this codebase
#   uses the natural Julia arr[row, col] convention. 
#   therefore TRANSPOSE
#   at the IO boundary so callers can always assume:
#       size(raster.data) == (rows, cols)
#       raster.data[r, c] is the value at row r, column c
#   and write_raster transposes back to GDAL's convention before writing.
#
#  

using ArchGDAL
const AG = ArchGDAL

"""
    RasterData

Lightweight container: 2D Float32 matrix indexed [row, col] + geotransform
+ projection + nodata sentinel. Carrying this around avoids re-reading the
DEM repeatedly in Monte Carlo loops.

The `geotransform` is the standard GDAL 6-element affine:
    gt[1] = top-left X
    gt[2] = pixel width  (positive)
    gt[3] = row rotation (usually 0)
    gt[4] = top-left Y
    gt[5] = column rotation (usually 0)
    gt[6] = pixel height (negative for north-up rasters)
"""
struct RasterData
    data::Matrix{Float32}
    geotransform::Vector{Float64}
    projection::String
    nodata::Float32
end

"""
    read_raster(path) -> RasterData

Read a raster file and return a RasterData with data indexed [row, col].
"""
function read_raster(path::AbstractString)
    isfile(path) ||
        throw(ArgumentError("Raster not found: '$path'. Check the path in your config (data.dem_raster or friction.prebuilt_path)."))
    return AG.read(path) do ds
        band = AG.getband(ds, 1)
        # AG.read(band) returns arr[col, row]; transpose so we end up
        # with the Julia-native [row, col] layout used everywhere downstream.
        raw = AG.read(band)
        arr = permutedims(raw)
        gt = AG.getgeotransform(ds)
        proj = AG.getproj(ds)
        nd = something(AG.getnodatavalue(band), -9999.0)
        RasterData(Float32.(arr), gt, proj, Float32(nd))
    end
end

"""
    write_raster(path, raster::RasterData)
    write_raster(path, data::AbstractMatrix, ref::RasterData)

Write a Float32 GeoTIFF using the reference raster's grid/projection.
Expects the input matrix to be indexed [row, col] (Julia convention); the
function transposes it back to GDAL's [col, row] before writing so the
resulting file aligns correctly with the source DEM in any GIS.
"""
function write_raster(path::AbstractString, raster::RasterData)
    write_raster(path, raster.data, raster)
end

function write_raster(path::AbstractString, data::AbstractMatrix, ref::RasterData)
    arr = Float32.(data)
    rows, cols = size(arr)
    # GDAL expects arr[col, row]; transpose back at the IO boundary.
    gdal_arr = permutedims(arr)
    AG.create(path;
              driver = AG.getdriver("GTiff"),
              width = cols, height = rows, nbands = 1,
              dtype = Float32,
              options = ["COMPRESS=DEFLATE", "TILED=YES"]) do ds
        AG.setgeotransform!(ds, ref.geotransform)
        AG.setproj!(ds, ref.projection)
        band = AG.getband(ds, 1)
        AG.setnodatavalue!(band, Float64(ref.nodata))
        AG.write!(band, gdal_arr)
    end
    return path
end

"""
    pixel_resolution(raster) -> Float64

Return the (positive) cell size of the raster, assumed square.
"""
pixel_resolution(r::RasterData) = abs(r.geotransform[2])

"""
    world_to_pixel(raster, x, y) -> (row, col)

Convert world (x, y) coordinates to 1-based (row, col) indices into the
[row, col] data matrix. Acts as the exact inverse of `pixel_to_world`:
any (x, y) inside a cell maps back to that cell's (row, col).
`gt[6]` is typically negative for north-up rasters, so a larger y
(further north) maps to a smaller row index, as expected.
"""
function world_to_pixel(r::RasterData, x::Real, y::Real)
    gt = r.geotransform
    # 0-based fractional cell offset from the raster origin
    px = (x - gt[1]) / gt[2]
    py = (y - gt[4]) / gt[6]
    # floor() so that any point inside a cell maps to that cell, making
    # this the exact inverse of pixel_to_world's centre-of-cell forward map.
    return (Int(floor(py)) + 1, Int(floor(px)) + 1)
end

"""
    pixel_to_world(raster, row, col) -> (x, y)

Inverse of `world_to_pixel`. Returns the centre of the (row, col) cell
in world coordinates.
"""
function pixel_to_world(r::RasterData, row::Int, col::Int)
    gt = r.geotransform
    x = gt[1] + (col - 1) * gt[2] + 0.5 * gt[2]
    y = gt[4] + (row - 1) * gt[6] + 0.5 * gt[6]
    return (x, y)
end

"""
    read_points(path; name_field="name") -> Vector{NamedTuple}

Read a point layer (.gpkg, .shp, ...). Returns vector of (name, x, y)
tuples in the layer's native CRS. Extra attribute fields are not returned.
"""
function read_points(path::AbstractString; name_field::AbstractString="name")
    out = NamedTuple{(:name, :x, :y), Tuple{String,Float64,Float64}}[]
    AG.read(path) do ds
        layer = AG.getlayer(ds, 0)
        AG.resetreading!(layer)
        for feat in layer
            geom = AG.getgeom(feat, 0)
            x = AG.getx(geom, 0)
            y = AG.gety(geom, 0)
            name = try
                String(AG.getfield(feat, name_field))
            catch
                "feature_$(AG.getfid(feat))"
            end
            push!(out, (name=name, x=x, y=y))
        end
    end
    return out
end

"""
    read_points_full(path; name_field="name") -> Vector{NamedTuple}

Like `read_points` but also returns the `type` attribute when present
(used to filter passes vs sites in `pass_frequency.jl`).
"""
function read_points_full(path::AbstractString;
                          name_field::AbstractString="name",
                          type_field::AbstractString="type")
    out = NamedTuple{(:name, :x, :y, :type),
                     Tuple{String,Float64,Float64,String}}[]
    AG.read(path) do ds
        layer = AG.getlayer(ds, 0)
        AG.resetreading!(layer)
        for feat in layer
            geom = AG.getgeom(feat, 0)
            x = AG.getx(geom, 0)
            y = AG.gety(geom, 0)
            name = try
                String(AG.getfield(feat, name_field))
            catch
                "feature_$(AG.getfid(feat))"
            end
            ttype = try
                String(AG.getfield(feat, type_field))
            catch
                ""
            end
            push!(out, (name=name, x=x, y=y, type=ttype))
        end
    end
    return out
end

"""
    find_point(points, name) -> NamedTuple or throws
"""
function find_point(points, name::AbstractString)
    idx = findfirst(p -> p.name == name, points)
    idx === nothing &&
        throw(ArgumentError("Point named '$name' not found. Available: $(map(p->p.name, points))"))
    return points[idx]
end
