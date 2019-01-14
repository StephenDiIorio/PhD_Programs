include("Particle.jl")
include("Field.jl")
using StaticArrays

# export ParticleCloud, getpos, getmom, pushcloud!


struct ParticleCloud
    numpart::UInt64
    particles::Vector{Particle}

    function ParticleCloud(numpart::UInt64, pos_f::Function, mom_f::Function)
        particles = Vector{Particle}(undef, numpart)

        for i = 1:numpart
            @inbounds particles[i] = Particle(pos_f(i), mom_f(i))
        end

        new(numpart, particles)
    end
end


function getpos(cloud::ParticleCloud)::Vector{MVector{3, Float64}}
    to_ret = Vector{MVector{3, Float64}}(undef, cloud.numpart)

    for i = 1:cloud.numpart
        @inbounds to_ret[i] = cloud.particles[i].pos
    end

    return to_ret
end


function getmom(cloud::ParticleCloud)::Vector{MVector{3, Float64}}
    to_ret = Vector{MVector{3, Float64}}(undef, cloud.numpart)

    for i = 1:cloud.numpart
        @inbounds to_ret[i] = cloud.particles[i].mom
    end

    return to_ret
end


function pushcloud!(cloud::ParticleCloud,
                    nx::UInt64, ny::UInt64, nz::UInt64, dt::Float64,
                    inv_step_size::SVector{3, Float64},
                    efield::Field, bfield::Field)
    for i = 1:cloud.numpart
        @inbounds local_e_force = weightedforce(efield,
                                                cloud.particles[i].pos,
                                                nx, ny, nz,
                                                inv_step_size)
        @inbounds local_b_force = weightedforce(bfield,
                                                cloud.particles[i].pos,
                                                nx, ny, nz,
                                                inv_step_size)
        @inbounds particlepush!(cloud.particles[i],
                                local_e_force, local_b_force,
                                dt)
    end
    return
end
