#Lennard-Jones cluster. 
#Record the first two moments of the coupling distance
#as well as the title response, mobility, pressure and kinetic energy.
#
#sine shear such that sin has one period in [-5,5]
#


PIDS=()

L=2
dt=0.0001
T=200000
burn=100000
interval=10000
seed=8156733675


etas=(0.5 0.25 0.1 0.05 0.025 0.01 0.005 0.0025)
betas=(2 4)

echo Launched on $(hostname) at $(date)

for beta in ${betas[@]}; do
	for eta in ${etas[@]}; do
		julia lj_coupling_mobility_sin.jl sync $L $eta $beta $burn $T $dt $interval $seed&
		PIDS+=($!)
		julia lj_coupling_mobility_sin.jl sticky $L $eta $beta $burn $T $dt $interval $seed&
		PIDS+=($!)
	done
done

for i in ${PIDS[@]}; do 
	wait $i 
	echo done $i
done
echo finished at $(date)