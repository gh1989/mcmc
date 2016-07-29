#ifndef OPTIONS_H
#define OPTIONS_H

namespace MCMC
{

class Options
{
    public:
        Options( int argc, char *argv[] );
        Options(){};

        //Options(const Options& )               = default;
        //Options( Options&& )                   = default;
        //Options& operator=( const Options& )   = default;
        //Options& operator=( Options&& )        = default;

        ~Options(){};
        
        unsigned int rng_seed() { return _rng_seed; }
        unsigned int mcmc_trials() { return _mcmc_trials; }
        unsigned int path_length() { return _path_length; }
        unsigned int extra_data_ratio() { return _extra_data_ratio; }
        unsigned int parallel_paths() { return _parallel_paths; }
        unsigned int burn() { return _burn; }
        unsigned int parameters() { return _parameters; }
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


    private:
        unsigned int _path_length      = 16;
        unsigned int _mcmc_trials      = 10000;
        unsigned int _burn             = 0;
        unsigned int _extra_data_ratio = 1;
        unsigned int _parallel_paths   = 25;
        unsigned int _rng_seed         = 0;
        unsigned int _parameters       = 1;

        double _parameter_proposal_variance  = 0.1;
        double _observation_noise_variance   = 0.1;
        double _trajectory_path_delta        = 0.02;
        double _diffusion_coefficient        = 0.001;
};
}

#endif // End OPTIONS_H
