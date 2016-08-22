#l (diffusion coefficient)
./bin/LFLangevin -a LFMCMC/real_sigma -N 1000000 -K 15 -M 1 -P 1 -l -o 0.00001 -p 0.001  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/real_sigma -N 1000000 -K 15 -M 1 -P 2 -l -o 0.0001 -p 0.001  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/real_sigma -N 1000000 -K 15 -M 1 -P 5 -l -o 0.001 -p 0.001  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/real_sigma -N 1000000 -K 15 -M 1 -P 10 -l -o 0.01 -p 0.001 -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/real_sigma -N 1000000 -K 15 -M 1 -P 20 -l -o 0.1 -p 0.001 -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/real_sigma -N 1000000 -K 15 -M 1 -P 40 -l -o 1.0 -p 0.001 -c 0.05 -i 0
#M (extra data ratio)
#./bin/LFLangevin -N 100000 -K 20 -M 1 -P 10 -o 0.01 p 0.0001 -c 0.05 -i 0 -B 0