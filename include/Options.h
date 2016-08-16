#ifndef OPTIONS_H
#define OPTIONS_H

#include <fstream>
#include <getopt.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>

namespace MCMC
{

class Options
{

    public:
    
        // Called with options.
        Options( int argc, char *argv[] );
        
        // Default options.
        Options();
        
        /*
        Options(Options&&)=default;
        Options& operator=(Options&&)=default;
        Options(const Options&)=default;
        Options& operator=(const Options&)=default;   
        ~Options(){std::cout<<"Options Desctructing"<<std::endl;};
        */
        
        ~Options()
        {
            std::cout<<"~Options() called from "<< this << std::endl;
        }
        
        bool infer_diffusion_parameters() { return _infer_diffusion_parameters; }
        bool infer_drift_parameters() { return _infer_drift_parameters; }
        
        
        int rng_seed() const { return _rng_seed; }
        int mcmc_trials() const { return _mcmc_trials; }
        int path_length() const { return _path_length; }
        int extra_data_ratio() const { return _extra_data_ratio; }
        int parallel_paths() const { return _parallel_paths; }
        int burn() const { return _burn; }
        int cutoff() const { return _cutoff; }
        int number_particles() const {return _number_particles;}
        
        double log_real_sigma() const {return _log_real_sigma;}
        double log_start_sigma() const {return _log_start_sigma;}
        double parameter_proposal_sigma() const { return _parameter_proposal_sigma; }
        double observation_noise_sigma() const { return _observation_noise_sigma; }
        double trajectory_path_delta() const { return _trajectory_path_delta; }
        double parameter_proposal_diffusion_sigma() const {return _parameter_proposal_diffusion_sigma; }
        
        std::string output_subfolder() { return _output_subfolder; }
        
        void set_number_particles( int n){ _number_particles = n; }
        void set_mcmc_trials( int n ){ _mcmc_trials = n; }
        void set_extra_data_ratio( int n ){ _extra_data_ratio = n; }
        void set_path_length( int n ){ _path_length = n; }
        void set_burn( int n ){ _burn = n; }
        void set_rng_seed( int n ){ _rng_seed = n; }
        void set_parallel_paths( int n ) { _parallel_paths=n; }
        void set_parameter_proposal_sigma( double d ){ _parameter_proposal_sigma=d; }
        void set_observation_noise_sigma( double d ){ _observation_noise_sigma=d; }
        void set_trajectory_path_delta( double d ){ _trajectory_path_delta=d; }
        void set_log_real_sigma(double d){_log_real_sigma=d;}
        void default_values();
        void set_output_subfolder(std::string &s)
        {
            std::cout<<"set_output_subfolder "<< s << std::endl;
            _output_subfolder = s; 
        }
        
        void print_header( std::ofstream &file );
        void print_options( std::ostream &o );

    private:
        bool _infer_diffusion_parameters; 
        bool _infer_drift_parameters;
    
        int _path_length;
        int _mcmc_trials;
        int _burn;
        int _extra_data_ratio;
        int _parallel_paths;
        int _rng_seed;
        int _cutoff;
        int _number_particles;
        
        double _parameter_proposal_sigma;
        double _observation_noise_sigma;
        double _trajectory_path_delta;
        double _log_real_sigma;
        double _log_start_sigma;
        double _parameter_proposal_diffusion_sigma;
        
        std::string _output_subfolder;
};


}

#endif // End OPTIONS_H
