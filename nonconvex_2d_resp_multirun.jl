print(string("Started: Two dimensional Non-Convex example with ", coupling, " coupling\n",
        "etas = ", etas, ", \n",
        "beta = ", beta, ", burn = ", burn_in_time, ", T = ", T, ", dt = ", dt, ",\n",  
        "sampling interval = ", sampling_interval, ", batchs = ", batches, ", total runs = ", num_runs, 
        ", seed = ", seed, "\n\n"))
flush(stdout)
pmap((batch,eta) -> pmap_run(batch, eta), repeat(1:batches, length(etas)), repeat(etas, inner = batches))