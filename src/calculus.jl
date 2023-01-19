#=
    Copyright (c) 2023 Louis Yu Lu, MIT License
=#

export
    deriv,
    deriv2,
    integral

function deriv(vs::Vector{T}, ts::Vector{T})::Vector{T} where {T<:AbstractFloat}
    n = length(vs)
    d1 = zeros(T, n)
    d1[1] = (vs[2] - vs[1]) / (ts[2] - ts[1])
    for i in 2:n-1
        d1[i] = (vs[i+1] - vs[i-1]) / (ts[i+1] - ts[i-1])
    end
    d1[n] = (vs[n] - vs[n-1]) / (ts[n] - ts[n-1])
    d1
end

function deriv(vs::Vector{T})::Vector{T} where {T<:AbstractFloat}
    n = length(vs)
    d1 = zeros(T, n)
    d1[1] = (vs[2] - vs[1])
    for i in 2:n-1
        d1[i] = 0.5 * (vs[i+1] - vs[i-1])
    end
    d1[n] = vs[n] - vs[n-1]
    d1
end

function deriv2(vs::Vector{T}, ts::Vector{T})::Vector{T} where {T<:AbstractFloat}
    n = length(vs)
    d2 = zeros(T, n)
    for i in 2:n-1
        d2[i] = ((vs[i+1] - vs[i]) / (ts[i+1] - ts[i]) - (vs[i] - vs[i-1]) / (ts[i] - ts[i-1])) * 2.0 / (ts[i+1] - ts[i-1])
    end
    d2[1] = d2[2] # take following value
    d2[n] = d2[n-1] # take previous value
    d2
end

function deriv2(vs::Vector{T})::Vector{T} where {T<:AbstractFloat}
    n = length(vs)
    d2 = zeros(T, n)
    for i in 2:n-1
        d2[i] = vs[i+1] - 2.0 * vs[i] + vs[i-1]
    end
    d2[1] = d2[2] # take following value
    d2[n] = d2[n-1] # take previous value
    d2
end

function integral(vs::Vector{T}, ts::Vector{T})::Vector{T} where {T<:AbstractFloat}
    n = length(vs)
    s = zeros(T, n)
    for i in 2:n
        s[i] = s[i-1] + 0.5 * (vs[i-1] + vs[i]) * (ts[i] - ts[i-1])
    end
    s
end

function integral(vs::Vector{T})::Vector{T} where {T<:AbstractFloat}
    n = length(vs)
    s = zeros(T, n)
    for i in 2:n
        s[i] = s[i-1] + 0.5 * (vs[i-1] + vs[i])
    end
    s
end
