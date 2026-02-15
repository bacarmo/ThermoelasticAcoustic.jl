module ThermoelasticAcoustic

# Write your package code here.
export my_sum

"""
    my_sum(collection)

Compute the sum of elements in `collection` using type-stable iteration.

# Arguments
- `collection`: An indexable collection of numeric elements

# Returns
- Sum of all elements with type matching `eltype(collection)`

# Examples
```jldoctest
julia> my_sum([1, 2, 3, 4])
10

julia> my_sum(1:10)
55
```
"""
function my_sum(collection)
    accumulator = zero(eltype(collection))
    @inbounds for index in eachindex(collection)
        accumulator += collection[index]
    end
    return accumulator
end

end
