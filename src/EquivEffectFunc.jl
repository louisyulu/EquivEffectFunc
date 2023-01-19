#=
    Copyright (c) 2023 Louis Yu Lu, MIT License
=#

module EquivEffectFunc

using Dierckx

include("calculus.jl")
include("control_points.jl")

export eef_extrema,
    eef_partition

"""
    eef_extrema(xs, ys, p_order, depend_on)->(trend, diff)

Find the equivalent effect function with control points on extrema, 
    return the trend and the difference between the original `ys`` and `trend``.
# Arguments
- `T<:AbstractFloat`: generic type
- `xs::Vector{T}`: input x values.
- `ys::Vector{T}`: input y values.
- `p_order::Int=3`: polynomial order of `trend` function
- `depend_on::Symbol=:second_deriv`: extrema on :ys_val, :first_deriv or :second_deriv
- `trend::Vector{T}`: returned trend values
- `diff::Vector{T}`: returned diff values
"""
function eef_extrema(xs::Vector{T}, ys::Vector{T}; p_order::Int=3, depend_on::Symbol=:second_deriv)::Tuple{Vector{T},Vector{T}} where {T<:AbstractFloat}
    if p_order < 0 || p_order > 4
        throw(error("Invalid p_order ($p_order)"))
    end
    hs = if depend_on == :ys_val
        ys
    elseif depend_on == :first_deriv
        deriv(ys, xs)
    elseif depend_on == :second_deriv
        deriv2(ys, xs)
    else
        throw(error("invalid ($depend_on), must be among :ys_val, :first_deriv and :second_deriv"))
    end
    cp = ControlPoints(xs, ys)
    find_extrema!(cp, hs)
    adjust_ends!(cp)
    spl = Spline1D(cp.xE, cp.gE, k=p_order + 1) # integral on control points
    trend = derivative(spl, xs)
    diff = ys - trend
    trend, diff
end

"""
    eef_extrema(ys, p_order, depend_on)->(trend, diff)

Find the equivalent effect function with control points on extrema, x values spread on unit intervals, 
    return the trend and the difference between the original `ys`` and `trend``.
# Arguments
- `T<:AbstractFloat`: generic type
- `ys::Vector{T}`: input y values.
- `p_order::Int=3`: polynomial order of `trend` function
- `depend_on::Symbol=:second_deriv`: extrema on :ys_val, :first_deriv or :second_deriv
- `trend::Vector{T}`: returned trend values
- `diff::Vector{T}`: returned diff values
"""
function eef_extrema(ys::Vector{T}; p_order::Int=3, depend_on::Symbol=:second_deriv)::Tuple{Vector{T},Vector{T}} where {T<:AbstractFloat}
    if p_order < 0 || p_order > 4
        throw(error("Invalid p_order ($p_order)"))
    end
    hs = if depend_on == :ys_val
        ys
    elseif depend_on == :first_deriv
        deriv(ys)
    elseif depend_on == :second_deriv
        deriv2(ys)
    else
        throw(error("invalid ($depend_on), must be among :ys_val, :first_deriv and :second_deriv"))
    end
    cp = ControlPoints(ys)
    find_extrema!(cp, hs)
    adjust_ends!(cp)
    spl = Spline1D(cp.xE, cp.gE, k=p_order + 1) # integral on control points
    trend = derivative(spl, cp.xs)
    diff = ys - trend
    trend, diff
end

"""
    eef_partition(xs, ys, p_order, depend_on)->(trend, diff)

Find the equivalent effect function with control points on partition, 
    return the trend and the difference between the original `ys`` and `trend``.
# Arguments
- `T<:AbstractFloat`: generic type
- `xs::Vector{T}`: input x values.
- `ys::Vector{T}`: input y values.
- `ys::Vector{T}`: input y values.
- `split_levels::Int`: partition splitting levels
- `p_order::Int=3`: polynomial order of `trend` function
- `depend_on::Symbol=:abs_2nd_deriv`: extrema on :abs_ys_val, :abs_1st_deriv or :abs_2nd_deriv
- `trend::Vector{T}`: returned trend values
- `diff::Vector{T}`: returned diff values
"""
function eef_partition(xs::Vector{T}, ys::Vector{T}, split_levels::Int; p_order::Int=3, depend_on::Symbol=:abs_2nd_deriv)::Tuple{Vector{T},Vector{T}} where {T<:AbstractFloat}
    if p_order < 0 || p_order > 4
        throw(error("Invalid p_order ($p_order)"))
    end
    if split_levels < 1 || 2^split_levels + 1 > length(xs)
        throw(error("Invalid split_levels ($split_levels)"))
    end
    hs = if depend_on == :abs_ys_val
        abs.(ys)
    elseif depend_on == :abs_1st_deriv
        abs.(deriv(ys, xs))
    elseif depend_on == :abs_2nd_deriv
        abs.(deriv2(ys, xs))
    else
        throw(error("Invalid ($depend_on), must be among :abs_ys_val, :abs_1st_deriv and :abs_2nd_deriv"))
    end
    cp = ControlPoints(xs, ys)
    partition!(cp, hs, split_levels)
    adjust_ends!(cp)
    spl = Spline1D(cp.xE, cp.gE, k=p_order + 1) # integral on control points
    trend = derivative(spl, xs)
    diff = ys - trend
    trend, diff
end

"""
    eef_partition(ys, p_order, depend_on)->(trend, diff)

Find the equivalent effect function with control points on partition, x values spread on unit intervals, 
    return the trend and the difference between the original `ys`` and `trend``.
# Arguments
- `T<:AbstractFloat`: generic type
- `ys::Vector{T}`: input y values.
- `split_levels::Int`: partition splitting levels
- `p_order::Int=3`: polynomial order of `trend` function
- `depend_on::Symbol=:abs_2nd_deriv`: extrema on :abs_ys_val, :abs_1st_deriv or :abs_2nd_deriv
- `trend::Vector{T}`: returned trend values
- `diff::Vector{T}`: returned diff values
"""
function eef_partition(ys::Vector{T}, split_levels::Int; p_order::Int=3, depend_on::Symbol=:abs_2nd_deriv)::Tuple{Vector{T},Vector{T}} where {T<:AbstractFloat}
    if p_order < 0 || p_order > 4
        throw(error("Invalid p_order ($p_order)"))
    end
    if split_levels < 1 || 2^split_levels + 1 > length(ys)
        throw(error("Invalid split_levels ($split_levels)"))
    end
    hs = if depend_on == :abs_ys_val
        abs.(ys)
    elseif depend_on == :abs_1st_deriv
        abs.(deriv(ys))
    elseif depend_on == :abs_2nd_deriv
        abs.(deriv2(ys))
    else
        throw(error("Invalid ($depend_on), must be among :abs_ys_val, :abs_1st_deriv and :abs_2nd_deriv"))
    end
    cp = ControlPoints(ys)
    partition!(cp, hs, split_levels)
    adjust_ends!(cp)
    spl = Spline1D(cp.xE, cp.gE, k=p_order + 1) # integral on control points
    trend = derivative(spl, cp.xs)
    diff = ys - trend
    trend, diff
end

end # module EquivEffectFunc
