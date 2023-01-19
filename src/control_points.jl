#=
    Copyright (c) 2023 Louis Yu Lu, MIT License
=#

export
    ControlPoints,
    find_extrema!,
    adjust_ends!,
    partition!

mutable struct ControlPoints{T<:AbstractFloat}
    xs::Vector{T}
    ys::Vector{T}
    unit_step::Bool # x change 1 on each step
    gs::Vector{T} # integration of ys
    indices::Vector{Int} # indices of control points
    xE::Vector{T} # x values of control points
    yE::Vector{T} # y values of control points
    gE::Vector{T} # integration values of control points
end

function ControlPoints(xs::Vector{T}, ys::Vector{T})::ControlPoints{T} where {T<:AbstractFloat}
    ControlPoints(xs, ys, false, T[], Int[], T[], T[], T[])
end

function ControlPoints(ys::Vector{T})::ControlPoints{T} where {T<:AbstractFloat}
    ControlPoints(collect(T(1):1:length(ys)), ys, true, T[], Int[], T[], T[], T[])
end

"""
find_extrema!(cp, hs)->count

Find control points on min and max values of hs, return the control points count.
# Arguments
- `T<:AbstractFloat`: generic type
- `cp::ControlPoints{T}`: struct of control points
- `hs::Vector{T}`: values on which to determine the extrema, can be the original values, first or second derivatives.
  the first derivative gives the inflexing point on the original function, and generate best results
- `count::Int`: returned the count of control points
"""
function find_extrema!(cp::ControlPoints{T}, hs::Vector{T})::Int where {T<:AbstractFloat}
    INDEX_NOT_SET = 0
    inc = INDEX_NOT_SET
    dec = INDEX_NOT_SET
    push!(cp.indices, 1)
    push!(cp.xE, cp.xs[1])
    push!(cp.yE, cp.ys[1])
    n = length(cp.ys)
    for i in 2:n
        hp = hs[i-1]
        hc = hs[i]
        if hp == hc
            continue
        elseif hp < hc
            # increasing
            if dec != INDEX_NOT_SET
                if i - 1 > dec
                    m = div(dec + i - 1, 2) # half way
                    push!(cp.indices, m)
                    push!(cp.xE, cp.xs[m])
                    push!(cp.yE, cp.ys[m])
                else
                    push!(cp.indices, dec)
                    push!(cp.xE, cp.xs[dec])
                    push!(cp.yE, cp.ys[dec])
                end
            end
            inc = i
            dec = INDEX_NOT_SET
        elseif hp > hc
            # decreasing
            if inc != INDEX_NOT_SET
                if i - 1 > inc
                    m = div(inc + i - 1, 2)
                    push!(cp.indices, m)
                    push!(cp.xE, cp.xs[m])
                    push!(cp.yE, cp.ys[m])
                else
                    push!(cp.indices, inc)
                    push!(cp.xE, cp.xs[inc])
                    push!(cp.yE, cp.ys[inc])
                end
            end
            dec = i
            inc = INDEX_NOT_SET
        end
    end
    push!(cp.indices, n)
    push!(cp.xE, cp.xs[n])
    push!(cp.yE, cp.ys[n])
    gs = if cp.unit_step
        integral(cp.ys)
    else
        integral(cp.ys, cp.xs)
    end
    for k in cp.indices
        push!(cp.gE, gs[k])
    end
    length(cp.yE)
end

function adjust_ends!(cp::ControlPoints{T}) where {T<:AbstractFloat}
    # even extension on both ends
    m = length(cp.yE)
    x = zeros(T, m + 2)
    y = zeros(T, m + 2)
    g = zeros(T, m + 2)
    for i in 1:m
        x[i+1] = cp.xE[i]
        y[i+1] = cp.yE[i]
        g[i+1] = cp.gE[i]
    end
    x[1] = 2 * cp.xE[1] - cp.xE[2]
    y[1] = cp.yE[2]
    g[1] = 2 * cp.gE[1] - cp.gE[2]

    x[m+2] = 2 * cp.xE[m] - cp.xE[m-1]
    y[m+2] = cp.yE[m-1]
    g[m+2] = 2 * cp.gE[m] - cp.gE[m-1]

    cp.xE = x
    cp.yE = y
    cp.gE = g
end

"""
partition!(cp, hs, split_levels)->count

Find control points on the density of hs, return the control points count.
# Arguments
- `T<:AbstractFloat`: generic type
- `cp::ControlPoints{T}`: struct of control points
- `hs::Vector{T}`: values on which to determine control points, can be the absolute values original function, of first or second derivatives.
  the absolute value of first derivative reflect the change rate, absolute value of second derivative reflect the change cuvature 
  which is suggested for better result, more control points are resulted in the dense area of hs
- `split_levels::Int`: levels of bi-partitions between x min and x max
- `count::Int`: returned the count of control points
"""
function partition!(cp::ControlPoints{T}, hs::Vector{T}, split_levels::Int)::Int where {T<:AbstractFloat}
    w = abs.(hs)  # absolute value as the weight
    mt = w .* cp.xs # moment
    ws = if cp.unit_step
        integral(w)
    else
        integral(w, cp.xs)
    end
    ms = if cp.unit_step
        integral(mt)
    else
        integral(mt, cp.xs)
    end
    m = length(cp.xs)
    pts = [1, 1, m]
    if ws[m] == 0
        t = (xs[m] - cp.xs[1]) / 2  #  0 weight, half way
    else
        t = ms[m] / ws[m]  # weight center
    end
    pts[2] = find_index(cp, t, 1, m)
    indices_list = Vector{Int}[pts] # partition indices at different level
    for split in 2:split_levels
        k = 2^split + 1
        pts = ones(Int, k)
        prev = indices_list[end]
        n = length(prev) - 1
        j = 2
        for i in 1:n
            dw = ws[prev[i+1]] - ws[prev[i]]
            if dw == T(0)
                t = (cp.xs[prev[i+1]] - cp.xs[prev[i]]) / 2 # 0 weight, half way
            else
                t = (ms[prev[i+1]] - ms[prev[i]]) / dw # weight center
            end
            pts[j] = find_index(cp, t, prev[i], prev[i+1])
            j += 1
            pts[j] = prev[i+1]
            j += 1
        end
        pts2 = Int[pts[1]]
        # remove duplicates
        for i in 2:length(pts)
            if pts[i] > pts[i-1]
                push!(pts2, pts[i])
            end
        end
        push!(indices_list, pts2)
    end
    cp.indices = indices_list[end]
    n = length(cp.indices)
    cp.xE = zeros(T, n)
    cp.yE = zeros(T, n)
    cp.gE = zeros(T, n)
    gs = if cp.unit_step
        integral(cp.ys)
    else
        integral(cp.ys, cp.xs)
    end
    for (i, k) in enumerate(cp.indices)
        cp.xE[i] = cp.xs[k]
        cp.yE[i] = cp.ys[k]
        cp.gE[i] = gs[k]
    end
    n
end

function find_index(cp::ControlPoints{T}, t::T, k1::Int, k2::Int)::Int where {T<:AbstractFloat}
    k = k1
    for i in k1:k2
        if cp.xs[i] <= t && t < cp.xs[i+1]
            k = i
            break
        end
    end
    k
end