#include "Options.h"

using namespace MCMC;

void Options::default_values()
{
    _path_length            = 16;
    _mcmc_trials            = 10000;
    _burn                   = 0;
    _extra_data_ratio       = 1;
    _parallel_paths         = 25;
    _rng_seed               = 0;
    _cutoff                 = 1;

    _parameter_proposal_sigma     = 0.05;
    _observation_noise_sigma      = 0.1;
    _trajectory_path_delta        = 0.001;
    _log_real_sigma               = log( 0.001 );
    _log_start_sigma              = 0.0; 
    _parameter_proposal_diffusion_sigma = 0.1;
    
    _infer_drift_parameters     = true;
    _infer_diffusion_parameters = false;
    _store_time_series          = true;
    _number_particles           = 100;
}

Options::Options()
{
    default_values();
}


Options::Options( int argc, char *argv[] )
{
    // The option.
    int opt;
    
    default_values();

    // Get options from command line.
    while( ( opt = getopt( argc, argv, ":R:K:P:N:B:M:c:o:p:g:l:d:D:i:Q:" ) ) != EOF ) 
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
            _parameter_proposal_sigma = atof(optarg);
            std::cout << "set _parameter_proposal_sigma = " << _parameter_proposal_sigma << std::endl;
            break;
            
            case 'o':
            _observation_noise_sigma = atof(optarg);
            std::cout << "set _observation_noise_sigma = " << _observation_noise_sigma << std::endl;
            break;
            
            case 'p':
            _trajectory_path_delta = atof(optarg);
            std::cout << "set _trajectory_path_delta = " << _trajectory_path_delta << std::endl;
            break;
            
            case 'g':
            _log_start_sigma = atof(optarg);
            std::cout << "set _log_start_sigma = "<< _log_start_sigma << std::endl;
            break;
            
            case 'l':
            _log_real_sigma = atof(optarg);
            std::cout<< "set _log_real_sigma ="<< _log_real_sigma << std::endl;
            break;
            
            case 'd':
            _parameter_proposal_diffusion_sigma = atof(optarg);
            std::cout << "set _parameter_proposal_diffusion_sigma = " << _parameter_proposal_diffusion_sigma << std::endl;
            break;
            
            case 'i':
            _infer_diffusion_parameters = atoi(optarg); 
            std::cout << "set _infer_diffusion_parameters = " << _infer_diffusion_parameters << std::endl;
            break;
            
            case 'D':
            _infer_drift_parameters = atoi(optarg);
            std::cout << "set _infer_drift_parameters = " << _infer_drift_parameters << std::endl;
            break;
            
            case 'Q':
            _number_particles = atoi(optarg);
            std::cout << "set _number_particles = " << _number_particles << std::endl;
            break;
            
        }
    }
}

void Options::print_header( std::ofstream &file )
{
    file << "#_log_real_sigma: " << _log_real_sigma <<std::endl;
    file << "#_observation_noise_sigma: " << _observation_noise_sigma <<std::endl;
    file << "#_path_length: "  << _path_length <<std::endl;
    file << "#_extra_data_ratio: " << _extra_data_ratio <<std::endl;
    file << "#_parallel_paths: "<< _parallel_paths <<std::endl;
    file << "#_parameter_proposal_sigma: " << _parameter_proposal_sigma <<std::endl;
}

void Options::print_options( std::ostream &o )
{
    o << "_output_subfolder = " << _output_subfolder <<std::endl;
    o << "_log_real_sigma = " << _log_real_sigma <<std::endl;
    o << "_observation_noise_sigma = " << _observation_noise_sigma <<std::endl;
    o << "_path_length = "  << _path_length <<std::endl;
    o << "_extra_data_ratio = " << _extra_data_ratio <<std::endl;
    o << "_parallel_paths = "<< _parallel_paths <<std::endl;
    o << "_parameter_proposal_sigma = " << _parameter_proposal_sigma <<std::endl;
    o << "_mcmc_trials = " << _mcmc_trials <<std::endl;
}
