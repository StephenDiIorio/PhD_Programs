include("Particle.jl")
include("ParticleCloud.jl")
include("Field.jl")
using StaticArrays
using Plots
using BenchmarkTools
pyplot()


function initpos(i::UInt64)::MVector{3, Float64}
    return @MVector fill(5, 3)
end


function initmom(i::UInt64)::MVector{3, Float64}
    return @MVector [1., 0., 0.]
end


function init_efield(nx::UInt64, ny::UInt64, nz::UInt64)::Array{Float64, 4}
    return zeros(Float64, 3, nx, ny, nz)
end


function init_bfield(nx::UInt64, ny::UInt64, nz::UInt64)::Array{Float64, 4}
    to_ret = zeros(Float64, 3, nx, ny, nz)
    to_ret[3, :, :, :] .= one(Float64)
    return to_ret
end


function main()
    t::Float64 = 0.0
    dt::Float64 = 0.1
    tmax::Float64 = 20. * pi
    npart::UInt64 = 1

    nx::UInt64 = 100
    ny::UInt64 = 100
    nz::UInt64 = 100

    stepsize = @SVector fill(0.1, 3)
    invstep = 1.0 ./ stepsize

    cloud = ParticleCloud(npart, initpos, initmom)

    efield = Field(nx, ny, nz, init_efield)
    bfield = Field(nx, ny, nz, init_bfield)

    c::UInt64 = 1
    data = Array{MVector{3, Float64}}(undef, round(Int64, tmax / dt, RoundUp))

    while t <= tmax
        pushcloud!(cloud, nx, ny, nz, dt, invstep, efield, bfield)

        data[c] = getpos(cloud)[1]

        c += 1
        t += dt
    end

    # println(data)
    plot(data[1], data[2], data[3])
    png("test")
end

main()

# dump(t)
# tune!(t)
# run(t)
