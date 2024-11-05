using Random
using LinearAlgebra 
using Statistics
using Distributions
using JLD
using Distributed

include("integrators.jl")
include("nonconvex_2d.jl")

#path for directory to save data and states to
data_dir = 
state_dir = 
params_path(L, batch, seed, eta, beta, burn_in_time, T, dt, sampling_interval, runs) = string("_L", L, "_batch", 
    batch, "_seed", seed, "_eta", eta, "_beta", beta, "_burn", burn_in_time, "_T", T, "_dt", dt, 
    "_interval", sampling_interval, "_numruns", runs, ".jld")


F(x) = [sin(x[2]), 0]


function R(x,y)
    return [x[1]*x[2] - y[1]*y[2], norm(x-y)]
end

R_output::Int64 = 2


U, grad_U = get_potential(L)

function drift(x, eta)
    return -grad_U(x) + eta*F(x)
end

runs_per_batch = num_runs ÷ batches

coupling_methods = Dict("indep"=>indep_noise, "sync"=>sync_noise, "reflect"=>reflect_noise, "sticky"=>sticky_noise)

function pmap_run(batch, eta)
    X_0 = [1., 0]

    Random.seed!(seed*batch)
    RNG = copy(Random.default_rng())

    N = floor(Int64, T/dt)
    states = ()
    R_traj = zeros(runs_per_batch, N ÷ sampling_interval, 2)
    R_avg = zeros(runs_per_batch, N ÷ sampling_interval, 2)

    print("η = $(eta), batch = $(batch) \n")
    flush(stdout)
    for j in 1:runs_per_batch
        R_traj[j, :, :], R_avg[j, :, :], state = compute_response(X_0, drift, coupling_methods[coupling], 
            R, R_output, eta, beta, burn_in_time, T, dt, sampling_interval, RNG)
        states = (states..., state)
    end
    save(data_dir*coupling*"_nonconvex2D_resp"*params_path(round(L; digits = 2), batch, seed, eta, beta, burn_in_time, 
        T, dt, sampling_interval, runs_per_batch), "R_traj", R_traj, "R_avg", R_avg)
    save(state_dir*coupling*"_nonconvex2D_resp_state"*params_path(round(L; digits = 2), batch, seed, eta, beta, 
        burn_in_time, T, dt, sampling_interval, runs_per_batch), "states", states)

    return nothing
end