#!/usr/bin/env julia
# Persistent, serial distance worker for fair_benchmark.py. Every request is
# original-path TAB simplified-path; stdout contains only one distance per row.
using FrechetDist
import FrechetDist.cg.point: npoint
import FrechetDist.cg.polygon: Polygon2F

function read_curve(path)
    lines = readlines(path)
    n = parse(Int, lines[1])
    n > 0 && length(lines) == n + 1 || error("Invalid curve: $path")
    curve = Polygon2F()
    for line in lines[2:end]
        xy = parse.(Float64, split(line))
        length(xy) == 2 && all(isfinite, xy) || error("Invalid point: $path")
        push!(curve, npoint(xy...))
    end
    return curve
end

for request in eachline(stdin)
    try
        paths = split(request, '\t')
        length(paths) == 2 || error("Expected two paths")
        distance = Float64(frechet_c_compute(read_curve(paths[1]),
                                           read_curve(paths[2])).leash)
        isfinite(distance) && distance >= 0 || error("Invalid distance")
        println(distance)
    catch err
        println("ERROR: ", sprint(showerror, err))
    end
    flush(stdout)
end
