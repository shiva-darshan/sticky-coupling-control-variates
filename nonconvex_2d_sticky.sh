#Multiple runs of sticky coupled dynamics in 2d non-convex potential
#500 runs over 5 batches for etas = 0.1, 0.05, 0.025, 0.01, 0.005
#


echo Launched on $(hostname) at $(date)
julia -p 25 -L /libre/darshans/scripts/sticky_nonconvex_2d_params.jl -L nonconvex_2d_resp_multirun_setup.jl nonconvex_2d_resp_multirun.jl
echo finished at $(date)