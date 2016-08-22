#o (observation noise)
./bin/LFLangevin -a LFMCMC/observation_noise_sigma -N 1000000 -K 15 -M 1 -P 10 -o 0.0005 -p 0.001 -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/observation_noise_sigma -N 1000000 -K 15 -M 1 -P 10 -o 0.001 -p 0.001  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/observation_noise_sigma -N 1000000 -K 15 -M 1 -P 10 -o 0.005 -p 0.001  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/observation_noise_sigma -N 1000000 -K 15 -M 1 -P 10 -o 0.01 -p 0.001  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/observation_noise_sigma -N 1000000 -K 15 -M 1 -P 10 -o 0.05 -p 0.001 -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/observation_noise_sigma -N 1000000 -K 15 -M 1 -P 10 -o 0.1 -p 0.001 -c 0.05 -i 0
#M (extra data ratio)
#./bin/LFLangevin -N 100000 -K 20 -M 1 -P 10 -o 0.01 p 0.0001 -c 0.05 -i 0 -B 0