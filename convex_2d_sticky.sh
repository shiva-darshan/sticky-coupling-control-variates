PID=()
beta=1
burn=100000
T=1000000
dt=0.005
interval=20000
num_runs=250
batchs=5
seeds=(9909724496 6924031646)

etas=(0.1 0.05 0.025 0.01 0.005)

for seed in ${seeds[@]}; do
	for eta in ${etas[@]}; do
		julia convex_2d_resp_multirun.jl sticky $eta $beta $burn $T $dt $interval $num_runs $batchs $seed &
		echo $eta $seed
	done
done
