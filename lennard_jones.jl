using LinearAlgebra

const e::Float64 = 1.
const s::Float64 = (2.)^(-1/6)
const a::Float64 = 1.
const M::Float64 = 5. #l^\infty distance at which coupling pot kicks in
const conf_pow::Float64 = 2.

v(r) = 4 * e *((s/r)^12 - (s/r)^6) #LJ
@inline @fastmath v_prime(r::Float64) = -4. * e * ((12. * s^12. / r^13.)- (6. * s^6. / r^7. ))
u_conf(q) = (max(0, abs(q[1]) - M)^conf_pow + max(0, abs(q[2]) - M)^conf_pow)/conf_pow #confining potential 
@fastmath grad_u_conf(q) = [max(0, abs(q[1]) - M)^(conf_pow-1) * sign(q[1]), max(0, abs(q[2]) - M)^(conf_pow-1) * sign(q[2])]

@fastmath grad_u_conf_comp(x)= max(0, abs(x) - M)^(conf_pow-1) * sign(x)

function U(q)
	d = length(q)
	n = floor(Int64, d/2)

	out = 0
	for i in 1:n
		for j in i:n
			if i == j
				out += a  * u_conf(q[2*i-1:2*i]) #confining potential
				out += v(norm(q[2*i-1:2*i])) #LJ potential
			else
				out += v(norm(q[2*i-1:2*i] - q[2*j-1:2*j]))
			end
		end
		
	end

	return out
end

function U2(q)
	d = length(q)
	n = floor(Int64, d/2)

	out = 0
	for i in 1:n
		out += a  * u_conf(q[2*i-1:2*i]) #confining potential
		for j in (i+1):n
			out += v(norm(q[2*i-1:2*i] - q[2*j-1:2*j]))
		end
		
	end

	return out
end

function grad_U!(gradient::Array{Float64,1}, q::Array{Float64,1})::Array{Float64,1}
	d = length(q)
	n = floor(Int64, d/2)
	leg = zeros(2)
	gradient .= 0.

	for i in 1:n
		for l = 1:2
			@inbounds leg[l] = q[2*i-2 + l]
		end
		r = norm(leg)
		for k in 1:2
			@inbounds gradient[2*i-2 + k] += a * grad_u_conf_comp(q[2*i-2 + k]) #confining force

			@inbounds gradient[2*i-2 + k] += leg[k] * v_prime(r)/r # interaction force
		end
		
		for j in (i+1):n
			for l = 1:2
				@inbounds leg[l] = q[2*i-2 + l] - q[2*j-2 + l]
			end
			r = norm(leg)
			for k in 1:2
				@fastmath @inbounds gradient[2*i-2 + k] += leg[k] * v_prime(r)/r # interaction force
				@fastmath @inbounds gradient[2*j-2 + k] -= leg[k] * v_prime(r)/r
			end
		end
	end

	return gradient
end

function grad_U2!(gradient::Array{Float64,1}, q::Array{Float64,1})::Array{Float64,1}
	#Potential without particle pinned to the origin
	#
	d = length(q)
	n = floor(Int64, d/2)
	leg = zeros(2)
	gradient .= 0.

	for i in 1:n
		for k in 1:2
			@inbounds gradient[2*i-2 + k] += a * grad_u_conf_comp(q[2*i-2 + k]) #confining force
		end
		
		for j in (i+1):n
			for l = 1:2
				@inbounds leg[l] = q[2*i-2 + l] - q[2*j-2 + l]
			end
			r = norm(leg)
			for k in 1:2
				@fastmath @inbounds gradient[2*i-2 + k] += leg[k] * v_prime(r)/r # interaction force
				@fastmath @inbounds gradient[2*j-2 + k] -= leg[k] * v_prime(r)/r
			end
		end
	end

	return gradient
end

function grad_U(q::Array{Float64,1})
	d = length(q)
	n = floor(Int64, d/2)
	leg = zeros(2)
	gradient = zeros(d)

	for i in 1:n
		for l = 1:2
			@inbounds leg[l] = q[2*i-2 + l]
		end
		r = norm(leg)
		for k in 1:2
			@inbounds gradient[2*i-2 + k] += a * grad_u_conf_comp(q[2*i-2 + k]) #confining force

			#particle pinned to the center
			@inbounds gradient[2*i-2 + k] += leg[k] * v_prime(r)/r # interaction force
		end
		
		for j in (i+1):n
			for l = 1:2
				@inbounds leg[l] = q[2*i-2 + l] - q[2*j-2 + l]
			end
			r = norm(leg)
			for k in 1:2
				@inbounds gradient[2*i-2 + k] += leg[k] * v_prime(r)/r # interaction force
				@inbounds gradient[2*j-2 + k] -= leg[k] * v_prime(r)/r
			end
		end
	end

	return gradient
end

function grad_U2(q::Array{Float64,1})
	#Potential without particle pinned to origin
	d = length(q)
	n = floor(Int64, d/2)
	leg = zeros(2)
	gradient = zeros(d)

	for i in 1:n
		for k in 1:2
			@inbounds gradient[2*i-2 + k] += a * grad_u_conf_comp(q[2*i-2 + k]) #confining force
		end
		
		for j in (i+1):n
			for l = 1:2
				@inbounds leg[l] = q[2*i-2 + l] - q[2*j-2 + l]
			end
			r = norm(leg)
			for k in 1:2
				@inbounds gradient[2*i-2 + k] += leg[k] * v_prime(r)/r # interaction force
				@inbounds gradient[2*j-2 + k] -= leg[k] * v_prime(r)/r
			end
		end
	end

	return gradient
end

function getX_0(L)
	e1 = [1, 0]
	e2 = [0.5, sqrt(3)/2]

	N = ceil(Int64, L)
	X_0 = []

	for i in -N:N, j in -N:N
	    if abs(i) + abs(j) == 0
	        continue 
	    end
	    
	    sample = i * e1 + j * e2
	    
	    if norm(sample) <= L
	        X_0 = [X_0; sample]
	    end
	end

	scaling = 0.1:0.001:3
	poss_init = zeros(length(scaling))

	for (i, scale) in zip(1:length(scaling), scaling)
	    poss_init[i] = U(scale * X_0)
	    
	end

	X_0 *= scaling[argmin(poss_init)]

	return X_0
end


function F_sin(q)
	d = length(q)
	shear_forcing = zeros(d)
	shear_forcing[1:2:end] = sin.(pi/M * q[2:2:end]) #was 2pi/M. Now particles with y-comp in [0,M] are pushed right and in [-M, 0] left

	return shear_forcing
end

function F_sin2(q)
	#Sin forcing normalized by system size
	#
	d = length(q)
	n = d ÷ 2
	shear_forcing = zeros(d)
	shear_forcing[1:2:end] = sin.(2pi/M * q[2:2:end])

	return shear_forcing ./ sqrt(n)
end

function F_sign(q)
	d = length(q)
	shear_forcing = zeros(d)
	shear_forcing[1:2:end] = sign.(q[2:2:end])

	return shear_forcing
end

const F_lin_norm::Float64 = 3.43036
function F_lin(q)
	d = length(q)
	forcing = zeros(d)
	for i = 1:2:d
		forcing[i] = q[i+1] / F_lin_norm
	end
	return forcing
end

function F_lin_unnormed(q)
	d = length(q)
	forcing = zeros(d)
	for i = 1:2:d
		forcing[i] = q[i+1]
	end
	return forcing
end

function comp_I(q, d)
	#from wikipedia entry on moment of inertia tensor
	I = zeros(2,2)

	I[1,1] = sum([q[2*k]^2 for k = 1:d])
	I[1,2] = sum([-q[2*k-1] * q[2k] for k = 1:d])
	I[2,1] = I[1,2]
	I[2,2] =  sum([q[2*k-1]^2 for k = 1:d])

	return I/d^2
end

function comp_I_off_diag(q,d)
	return sum([-q[2*k-1] * q[2k] for k = 1:d])
end

function max_dist(q,d)
	return maximum([norm(q[2*k-1:2k]) for k = 1:d])
end

function min_dist(q,d)
	return minimum([norm(q[2*k-1:2k]) for k = 1:d])
end

