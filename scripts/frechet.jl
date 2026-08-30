#!/usr/bin/env julia
"""
Fréchet distance wrapper — supports both single-path and batch modes.

Single-path mode (legacy):
  frechet.jl <id>                       Compare data/<id>/original.txt with
                                        data/<id>/simplify.txt.
  frechet.jl --in <id>                  Same as positional <id>.
  frechet.jl --in <id> --path <file>    Compare original against an explicit
                                        simplified-curve file.
  frechet.jl --raw                      Print only the raw distance number.

Batch mode (benchmark speedup):
  frechet.jl --id <id> --batch <path1> --batch <path2> ... --batch <pathN>
                                        Reads original once, computes Frechet against
                                        all N batch paths in a single Julia session,
                                        then prints one line per path:
                                          <basename1>: <distance1>
                                          <basename2>: <distance2>
                                          ...
                                        All in ~1s instead of N × 7s overhead.

Usage:
  julia frechet.jl --id 42 --batch data/42/DP.txt --batch data/42/simplify_against_DP_0.7_123.4.txt --batch data/42/simplify_against_DOTS_0.8_456.7.txt --raw
"""

using FrechetDist
using ArgParse

import FrechetDist.cg.point: npoint
import FrechetDist.cg.polygon: Polygon2F

function read_curve(filename)
    lines = readlines(filename)
    isempty(lines) && error("Empty file: $filename")

    n = try
        parse(Int, strip(lines[1]))
    catch
        error("Expected first line to be integer N in $filename")
    end

    if length(lines) < n + 1
        error("Expected $n points but only $(length(lines) - 1) lines")
    end

    poly = Polygon2F()
    for i in 1:n
        parts = split(strip(lines[i + 1]))
        if length(parts) < 2
            error("Malformed point at line $(i+1)")
        end
        x = parse(Float64, parts[1])
        y = parse(Float64, parts[2])
        push!(poly, npoint(x, y))
    end

    return poly
end

function run_batch(id::Int, batch_paths::AbstractVector, raw::Bool)
    repo_root = dirname(dirname(@__FILE__))
    original_path = joinpath(repo_root, "data", string(id), "original.txt")

    if !isfile(original_path)
        println(stderr, "Original file not found: $original_path")
        return 1
    end

    P = read_curve(original_path)

    for path in batch_paths
        if !isfile(path)
            println(stderr, "MISSING: $path")
            continue
        end
        try
            Q = read_curve(path)
            mw = frechet_c_compute(P, Q)
            dist = Float64(mw.leash)
            println("$(basename(path)): $dist")
            flush(stdout)
        catch e
            println(stderr, "ERROR ($path): $e")
        end
    end

    return 0
end

function run_single(id::Int, path_arg::String, raw::Bool)
    repo_root = dirname(dirname(@__FILE__))
    base = joinpath(repo_root, "data", string(id))
    original_path = joinpath(base, "original.txt")

    if !isfile(original_path)
        println(stderr, "Original file not found: $original_path")
        return 1
    end

    targets = if !isempty(path_arg)
        [path_arg]
    else
        [joinpath(base, "simplify.txt")]
    end

    P = read_curve(original_path)

    for p in targets
        if !isfile(p)
            println("  MISSING: $p")
            continue
        end
        try
            Q = read_curve(p)
            mw = frechet_c_compute(P, Q)
            dist = mw.leash
            if raw
                println(dist)
            else
                println("Frechet distance: $dist")
            end
            flush(stdout)
        catch e
            println(stderr, "  ERROR ($p): $e")
        end
    end

    return 0
end

function main()
    s = ArgParseSettings()
    s.description = "Compute continuous Fréchet distances."

    @add_arg_table s begin
        "id"
            help = "trajectory id (single-path legacy mode)"
            arg_type = Int
            required = false
        "--in"
            help = "trajectory id (alternative to positional)"
            arg_type = Int
        "--path", "-p"
            help = "path to simplified curve file (single-path mode)"
            arg_type = String
            default = ""
        "--raw"
            help = "raw output only"
            action = :store_true
        "--id"
            help = "trajectory id (batch mode)"
            arg_type = Int
            default = 0
        "--batch", "-b"
            help = "path to a simplified curve (batch mode, repeatable)"
            arg_type = String
            action = :append_arg
    end

    parsed = parse_args(s)

    # Batch mode: --id is set with --batch
    batch_id = parsed["id"]
    batch_paths_list = parsed["batch"]
    has_batch_id = batch_id != 0 && !isempty(batch_paths_list)

    if has_batch_id
        return run_batch(batch_id, batch_paths_list, parsed["raw"])
    else
        # Legacy single-path mode: use positional "id" or --in, NOT --id
        id = parsed["id"] !== nothing && parsed["id"] != 0 ? parsed["id"] : parsed["in"]
        if id === nothing
            println(stderr, "Error: trajectory id is required")
            println(stderr, "  Batch mode: --id <id> --batch <path1> --batch <path2> ...")
            println(stderr, "  Single mode: <id> [--path <file>] [--raw]")
            return 1
        end
        return run_single(id, parsed["path"], parsed["raw"])
    end
end

exit(main())
