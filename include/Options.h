#ifndef OPTIONS_H
#define OPTIONS_H

#include <getopt.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>

namespace MCMC
{

template <class Dynamics>
class Options
{
    typedef typename Dynamics::ParameterType ParameterType;

    public:
        Options( int argc, char *argv[] );

        Options()
        {
            _parameters = Dynamics::default_parameters(*this);
        }

        ~Options(){}
        
        unsigned int rng_seed() { return _rng_seed; }
        unsigned int mcmc_trials() { return _mcmc_trials; }
        unsigned int path_length() { return _path_length; }
        unsigned int extra_data_ratio() { return _extra_data_ratio; }
        unsigned int parallel_paths() { return _parallel_paths; }
        unsigned int burn() { return _burn; }
        unsigned int cutoff() { return _cutoff; }
        double parameter_proposal_variance() { return _parameter_proposal_variance; }
        double observation_noise_variance() { return _observation_noise_variance; }
        double trajectory_path_delta() { return _trajectory_path_delta; }
        double diffusion_coefficient() { return _diffusion_coefficient; }
    

        void set_mcmc_trials( int n ){ _mcmc_trials = n; }
        void set_extra_data_ratio( int n ){ _extra_data_ratio = n; }
        void set_path_length( int n ){ _path_length = n; }
        void set_burn( int n ){ _burn = n; }
        void set_rng_seed( int n ){ _rng_seed = n; }
        void set_parallel_paths( int n ) { _parallel_paths=n; }
        void set_parameter_proposal_variance( double d ){ _parameter_proposal_variance=d; }
        void set_observation_noise_variance( double d ){ _observation_noise_variance=d; }
        void set_trajectory_path_delta( double d ){ _trajectory_path_delta=d; }
        void set_diffusion_coefficient( double d ){ _diffusion_coefficient=d; }

        void set_parameters( ParameterType p ){ _parameters = p; }
        ParameterType parameters() { return _parameters; }

    private:
        unsigned int _path_length      = 16;
        unsigned int _mcmc_trials      = 10000;
        unsigned int _burn             = 0;
        unsigned int _extra_data_ratio = 1;
        unsigned int _parallel_paths   = 25;
        unsigned int _rng_seed         = 0;
        unsigned int _cutoff           = 1;

        double _parameter_proposal_variance  = 0.1;
        double _observation_noise_variance   = 0.1;
        double _trajectory_path_delta        = 0.02;
        double _diffusion_coefficient        = 0.001;

        ParameterType _parameters;

};

}

using namespace MCMC;

template <class Dynamics>
Options<Dynamics>::Options( int argc, char *argv[] )
{
    // The option.
    int opt;
    
    // Get options from command line.
    while( ( opt = getopt( argc, argv, ":R:K:P:N:B:M:c:o:p:d:" ) ) != EOF ) 
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
            _diffusion_coefficient = atof(optarg);
            std::cout << "set _diffusion_coefficient = " << _diffusion_coefficient << std::endl;
            break;

        }
    }
    std::cout<<"Here 1."<<std::endl;

    _parameters = Dynamics::default_parameters(*this);

    std::cout<<"Here 2."<<std::endl;
}

#endif // End OPTIONS_H
