#M (extra data ratio)
./bin/LFLangevin -a LFMCMC/extra_data_ratio -N 1000000 -K 15 -M 1 -P 10 -o 0.01 -p 0.001  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/extra_data_ratio -N 1000000 -K 15 -M 2 -P 10 -o 0.01 -p 0.0005  -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/extra_data_ratio -N 1000000 -K 15 -M 3 -P 10 -o 0.01 -p 0.00033333 -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/extra_data_ratio -N 1000000 -K 15 -M 4 -P 10 -o 0.01 -p 0.00025 -c 0.05 -i 0 
./bin/LFLangevin -a LFMCMC/extra_data_ratio -N 1000000 -K 15 -M 5 -P 10 -o 0.01 -p 0.0002 -c 0.05 -i 0
./bin/LFLangevin -a LFMCMC/extra_data_ratio -N 1000000 -K 15 -M 6 -P 10 -o 0.01 -p 0.00016666 -c 0.05 -i 0