mutable struct MutableFloat
    value::Float64
end

function set_all!(
	vec::Vector, value::T = 0
) where {T <: Real}
	for ii = 1:length(vec)
		if typeof(vec[ii]) <: Vector
			vec[ii] = set_all!(vec[ii], value)
		elseif typeof(vec[ii]) <: T
			vec[ii] = value
		end
	end
    return vec
end

# Assumes dst and src are vectors of vectors with the same 
# structure/length/types
function recursive_copy!(
	dst::Vector, 
	src::Vector
)
	@assert length(dst) == length(src)
	for ii = 1:length(dst)
		if typeof(dst[ii]) <: Vector
			@assert typeof(dst[ii]) == typeof(src[ii])
			recursive_copy!(dst[ii], src[ii])
		else
			dst[ii] = src[ii]
		end
	end
end