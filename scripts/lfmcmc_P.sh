#P (path length)
./bin/LFLangevin -a LFMCMC/path_length -N 1000000 -K 15 -M 1 -P 1 -o 0.01 -p 0.001  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/path_length -N 1000000 -K 15 -M 1 -P 2 -o 0.01 -p 0.001  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/path_length -N 1000000 -K 15 -M 1 -P 5 -o 0.01 -p 0.001  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/path_length -N 1000000 -K 15 -M 1 -P 10 -o 0.01 -p 0.001 -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/path_length -N 1000000 -K 15 -M 1 -P 20 -o 0.01 -p 0.001 -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/path_length -N 1000000 -K 15 -M 1 -P 40 -o 0.01 -p 0.001 -c 0.05 -i 0
#M (extra data ratio)
#./bin/LFLangevin -N 100000 -K 20 -M 1 -P 10 -o 0.01 p 0.0001 -c 0.05 -i 0 -B 0