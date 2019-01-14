using StaticArrays

# export Field, weightedforce


struct Field
    field::Array{Float64, 4}

    function Field(nx::UInt64, ny::UInt64, nz::UInt64, f::Function)
        new(f(nx, ny, nz))
    end
end


function weightedforce(field::Field, pos::MVector{3, Float64},
                       nx::UInt64, ny::UInt64, nz::UInt64,
                       inv_step_size::SVector{3, Float64})::MVector{3, Float64}
    # This algorithm was written dependent on a language that uses a 0-base
    # index. Julia is a 1-based index system. As such, 1 needs to be added to
    # ii in order for its index to be in bounds for a Julia array.
    qoverm = one(Float64)

    localforce = zeros(MVector{3, Float64})
    volumes = zeros(MArray{Tuple{2, 2, 2}, Float64})

    weight = copy(pos)
    weight .*= inv_step_size
    ii = map(x -> round(Int64, x, RoundDown), weight) # Add one when necessary
    weight .-= map(x -> convert(Float64, x), ii)

    for i = 1:size(volumes)[1], j = 1:size(volumes)[2], k = 1:size(volumes)[3]
        # Subtract 1 from i, j, k so they work with the calculation.
        @inbounds volumes[i, j, k] = abs(((one(Float64) - (i - 1)) - weight[1]) * ((one(Float64) - (j - 1)) - weight[2]) * ((one(Float64) - (k - 1)) - weight[3]))
    end

    for i = 1:size(volumes)[1], j = 1:size(volumes)[2], k = 1:size(volumes)[3]
        # Adding 1 to ii and subtracting 1 from i, j, k cancel out here
        i0 = ii[1] + i
        i1 = ii[2] + j
        i2 = ii[3] + k
        # Adjust the index value checks to match Julia array indexing
        if i0 <= nx && i1 <= ny && i2 <= nz &&
           i0 > 0 && i1 > 0 && i2 > 0
            @. @inbounds localforce += volumes[i, j, k] * field.field[:, i0, i1, i2]
        else
            println("WARNING: Evaluating field outside specified domain. (It treats this as if there is no field at this position.)")
        end
    end

    return localforce .* qoverm
end
