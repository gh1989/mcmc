#include "Options.h"
#include "lest/lest.hpp"

using namespace std;
using namespace MCMC;

Options opts;

const lest::test specification[] =
{

    CASE( "Options_rng_seed" )
    {
       opts.set_rng_seed(42);
       EXPECT( 42 == opts.rng_seed() );
    },
    
    CASE( "Options_mcmc_trials" )
    {
        opts.set_mcmc_trials(144758);
        EXPECT( 144758 == opts.mcmc_trials() );
    },

    CASE( "Options_path_length" )
    {
        opts.set_path_length(11234);
        EXPECT( 11234 == opts.path_length() ); 

    },

    CASE( "Options_extra_data_ratio" )
    {
        opts.set_extra_data_ratio(42);
        EXPECT( 42 == opts.extra_data_ratio() ); 

    },

    CASE( "Options_parallel_paths" )
    {
        opts.set_parallel_paths(42);
        EXPECT( 42 == opts.parallel_paths() ); 

    },

    CASE( "Options_burn" )
    {
        opts.set_burn(10001);
        EXPECT( 10001 == opts.burn() ); 
    },


    CASE( "Options_parameter_proposal_variance" )
    {
        opts.set_parameter_proposal_sigma(50.50);
        EXPECT( 50.50 == opts.parameter_proposal_sigma() ); 
    },

    CASE( "Options_observation_noise_variance" )
    {
        opts.set_observation_noise_sigma(0.4111);
        EXPECT( 0.4111 == opts.observation_noise_sigma() ); 
    },

    CASE( "Options_trajectory_path_delta" )
    {
        opts.set_trajectory_path_delta(0.42);
        EXPECT( 0.42 == opts.trajectory_path_delta() ); 
    },


};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv /*, std::cout */ );
}
