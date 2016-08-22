#K (parallel paths)
#nice -n 15 ./bin/LFLangevin -a LFMCMC/parallel_paths -N 1000000 -K 1 -M 1 -P 10 -o 0.01 -p 0.001  -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/parallel_paths -N 1000000 -K 2 -M 1 -P 10 -o 0.01 -p 0.001  -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/parallel_paths -N 1000000 -K 5 -M 1 -P 10 -o 0.01 -p 0.001  -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/parallel_paths -N 1000000 -K 10 -M 1 -P 10 -o 0.01 -p 0.001 -c 0.05 -i 0 
nice -n 15 ./bin/LFLangevin -a LFMCMC/parallel_paths -N 1000000 -K 20 -M 1 -P 10 -o 0.01 -p 0.001 -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/parallel_paths -N 1000000 -K 25 -M 1 -P 10 -o 0.01 -p 0.001 -c 0.05 -i 0
#l (diffusion coefficient)
#nice -n 15 ./bin/LFLangevin -a LFMCMC/real_sigma -N 1000000 -K 15 -M 1 -P 1 -l -o 0.00001 -p 0.001  -c 0.05 -i 0
#nice -n 15 ./bin/LFLangevin -a LFMCMC/real_sigma -N 1000000 -K 15 -M 1 -P 2 -l -o 0.0001 -p 0.001  -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/real_sigma -N 1000000 -K 15 -M 1 -P 5 -l -o 0.001 -p 0.001  -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/real_sigma -N 1000000 -K 15 -M 1 -P 10 -l -o 0.01 -p 0.001 -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/real_sigma -N 1000000 -K 15 -M 1 -P 20 -l -o 0.1 -p 0.001 -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/real_sigma -N 1000000 -K 15 -M 1 -P 40 -l -o 1.0 -p 0.001 -c 0.05 -i 0
#P (path length)
#nice -n 15 ./bin/LFLangevin -a LFMCMC/path_length -N 1000000 -K 15 -M 1 -P 1 -o 0.01 -p 0.001  -c 0.05 -i 0
#nice -n 15 ./bin/LFLangevin -a LFMCMC/path_length -N 1000000 -K 15 -M 1 -P 2 -o 0.01 -p 0.001  -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/path_length -N 1000000 -K 15 -M 1 -P 5 -o 0.01 -p 0.001  -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/path_length -N 1000000 -K 15 -M 1 -P 10 -o 0.01 -p 0.001 -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/path_length -N 1000000 -K 15 -M 1 -P 20 -o 0.01 -p 0.001 -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/path_length -N 1000000 -K 15 -M 1 -P 40 -o 0.01 -p 0.001 -c 0.05 -i 0
#o (observation noise)
nice -n 15 ./bin/LFLangevin -a LFMCMC/observation_noise_sigma -N 1000000 -K 15 -M 1 -P 10 -o 0.0005 -p 0.001 -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/observation_noise_sigma -N 1000000 -K 15 -M 1 -P 10 -o 0.001 -p 0.001  -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/observation_noise_sigma -N 1000000 -K 15 -M 1 -P 10 -o 0.005 -p 0.001  -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/observation_noise_sigma -N 1000000 -K 15 -M 1 -P 10 -o 0.01 -p 0.001  -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/observation_noise_sigma -N 1000000 -K 15 -M 1 -P 10 -o 0.05 -p 0.001 -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/observation_noise_sigma -N 1000000 -K 15 -M 1 -P 10 -o 0.1 -p 0.001 -c 0.05 -i 0
#M (extra data ratio)
nice -n 15 ./bin/LFLangevin -a LFMCMC/extra_data_ratio -N 1000000 -K 15 -M 1 -P 10 -o 0.01 -p 0.001  -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/extra_data_ratio -N 1000000 -K 15 -M 2 -P 10 -o 0.01 -p 0.0005  -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/extra_data_ratio -N 1000000 -K 15 -M 3 -P 10 -o 0.01 -p 0.00033333 -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/extra_data_ratio -N 1000000 -K 15 -M 4 -P 10 -o 0.01 -p 0.00025 -c 0.05 -i 0 
nice -n 15 ./bin/LFLangevin -a LFMCMC/extra_data_ratio -N 1000000 -K 15 -M 5 -P 10 -o 0.01 -p 0.0002 -c 0.05 -i 0
nice -n 15 ./bin/LFLangevin -a LFMCMC/extra_data_ratio -N 1000000 -K 15 -M 6 -P 10 -o 0.01 -p 0.00016666 -c 0.05 -i 0