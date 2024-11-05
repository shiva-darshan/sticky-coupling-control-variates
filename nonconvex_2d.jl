function get_potential(L)
	"""
	Generate locally non-convex potential with quadratic confinement outside
	an 2L by 2L box.
	"""

	a::Float64 = 2pi/L
	@inline function U(x)
		if norm(x, Inf) <= L
        	return (1 - cos(a*x[1]))*(1 - cos(a*x[2]))
	    else
	    	return a^2 * sum(max.(0, abs.(x) .-L).^2)/2
	    end
	end

	@inline function grad_U(x)
		if norm(x, Inf) <= L
        	return a .* [sin(a*x[1])*(1 - cos(a*x[2])), (1 - cos(a*x[1]))*sin(a*x[2])]
	    elseif abs(x[1]) > L
	        if abs(x[2]) > L
	            return  a^2 * (x - L*sign.(x))
	        else
	            return a^2 * [x[1] - L*sign(x[1]), 0]
	        end
	    elseif abs(x[2]) > L
	        return  a^2 * [0, x[2] - L*sign(x[2])]
	    end
	end

	return U, grad_U
end