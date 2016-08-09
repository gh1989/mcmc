#include "Options.h"

using namespace MCMC;

Options::Options( int argc, char *argv[] )
{
    // The option.
    int opt;
    
    // Get options from command line.
    while( ( opt = getopt( argc, argv, ":R:K:P:N:B:M:c:o:p:d:z:" ) ) != EOF ) 
    {
    switch (opt)
        {
            case 'R':
            _rng_seed = atoi(optarg);
            std::cout << "set _rng_seed = " << _rng_seed << std::endl;
            break;
            
            case 'K':
            _parallel_paths = atoi(optarg);
            std::cout << "set _parallel_paths = " << _parallel_paths << std::endl;
            break;
            
            case 'P':
            _path_length = atoi(optarg);
            std::cout << "set _path_length = " << _path_length << std::endl;
            break;
            
            case 'N':
            _mcmc_trials = atoi(optarg);
            std::cout << "set _mcmc_trials = " << _mcmc_trials << std::endl;
            break;
            
            case 'B':
            _burn = atoi(optarg);
            std::cout << "set _burn = " << _burn << std::endl;
            break;

            case 'M':
            _extra_data_ratio = atoi(optarg);
            std::cout<< "set _extra_data_ratio = " << _extra_data_ratio << std::endl;
            break;

            case 'c':
            _parameter_proposal_variance = atof(optarg);
            std::cout << "set _parameter_proposal_variance = " << _parameter_proposal_variance << std::endl;
            break;
            
            case 'o':
            _observation_noise_variance = atof(optarg);
            std::cout << "set _observation_noise_variance = " << _observation_noise_variance << std::endl;
            break;
            
            case 'p':
            _trajectory_path_delta = atof(optarg);
            std::cout << "set _trajectory_path_delta = " << _trajectory_path_delta << std::endl;
            break;
            
            case 'd':
            _variance_sigma = atof(optarg);
            std::cout << "set _variance_sigma = " << _variance_sigma << std::endl;
            break;
            
            case 'z':
            _zeta_sigma = atof(optarg);
            std::cout<< "set _zeta_sigma = " << _zeta_sigma << std::endl;
            break;

        }
    }

}