#K (parallel paths)
./bin/LFLangevin -a LFMCMC/parallel_paths -N 1000000 -K 1 -M 1 -P 10 -o 0.01 -p 0.001  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/parallel_paths -N 1000000 -K 2 -M 1 -P 10 -o 0.01 -p 0.001  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/parallel_paths -N 1000000 -K 5 -M 1 -P 10 -o 0.01 -p 0.001  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/parallel_paths -N 1000000 -K 10 -M 1 -P 10 -o 0.01 -p 0.001 -c 0.05 -i 0 
./bin/LFLangevin -a LFMCMC/parallel_paths -N 1000000 -K 20 -M 1 -P 10 -o 0.01 -p 0.001 -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/parallel_paths -N 1000000 -K 25 -M 1 -P 10 -o 0.01 -p 0.001 -c 0.05 -i 0