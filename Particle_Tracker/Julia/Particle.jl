using LinearAlgebra
using StaticArrays

# export Particle, particlepush!


mutable struct Particle
    pos::MVector{3, Float64}
    mom::MVector{3, Float64}
end


function particlepush!(part::Particle,
                       eforce::MVector{3, Float64},
                       bforce::MVector{3, Float64},
                       dt::Float64)
    pushmomentum!(part, eforce, dt / 0.5)
    lorentz!(part, bforce, dt)
    pushmomentum!(part, eforce, dt / 0.5)
    pushpos!(part, dt)
    return
end


function pushmomentum!(part::Particle,
                       local_e_force::MVector{3, Float64},
                       dt::Float64)
    @. part.mom += local_e_force * dt
    return
end


function pushpos!(part::Particle, dt::Float64)
    mom2 = dot(part.mom, part.mom)
    gamma = 1.0 / sqrt(1.0 + mom2)
    @. part.pos += part.mom * dt * gamma
    return
end


function lorentz!(part::Particle,
                  local_b_force::MVector{3, Float64},
                  dt::Float64)
    loc_b2 = dot(local_b_force, local_b_force)

    if loc_b2 != 0.0
        mom2 = dot(part.mom, part.mom)
        gamma = 1.0 / sqrt(1.0 + mom2)

        # Pre allocate arrays
        tt = MVector{3}(Vector{Float64}(undef, 3))
        ss = MVector{3}(Vector{Float64}(undef, 3))
        vperp = MVector{3}(Vector{Float64}(undef, 3))
        vstar = MVector{3}(Vector{Float64}(undef, 3))

        @. tt = local_b_force * dt * 0.5
        tt2 = dot(tt, tt)

        @. ss = 2.0 * tt / (1.0 + tt2)

        @. vperp = part.mom * (1.0 - local_b_force / loc_b2)
        c = cross(vperp, tt)
        @. vstar = vperp + c * gamma

        c = cross(vstar, ss)
        @. part.mom += c * gamma
    end
    return
end
