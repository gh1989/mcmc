#include "Options.h"
#include "OUPathScheme.h"
#include "lest/lest.hpp"

using namespace std;
using namespace MCMC;

const lest::test specification[] =
{
    CASE( "Options_Options" )
    {
        int argc=0;
        char **argv;
        Options opts(argc, argv);

    },

    CASE( "Options_rng_seed" )
    {
        int argc=3;
        char *argv[3];
        argv[1] = (char*)"-R";
        argv[2] = (char*)"42";
        Options opts(argc, argv);

        EXPECT( 42 == opts.rng_seed() );
    },

    CASE( "Options_mcmc_trials" )
    {
        int argc=3;
        char *argv[3];
        argv[1] = (char*)"-N";
        argv[2] = (char*)"144758";
        Options opts(argc, argv);

        EXPECT( 144758 == opts.mcmc_trials() );
    },

    CASE( "Options_path_length" )
    {
        int argc=3;
        char *argv[3];
        argv[1] = (char*)"-P";
        argv[2] = (char*)"11234";
        Options opts(argc, argv);

        EXPECT( 11234 == opts.path_length() ); 

    },

    CASE( "Options_extra_data_ratio" )
    {
        int argc=3;
        char *argv[3];
        argv[1] = (char*)"-M";
        argv[2] = (char*)"42";
        Options opts(argc, argv);

        EXPECT( 42 == opts.extra_data_ratio() ); 

    },

    CASE( "Options_parallel_paths" )
    {
        int argc=3;
        char *argv[3];
        argv[1] = (char*)"-K";
        argv[2] = (char*)"42";
        Options opts(argc, argv);

        EXPECT( 42 == opts.parallel_paths() ); 

    },

    CASE( "Options_burn" )
    {
        int argc=3;
        char *argv[3];
        argv[1] = (char*)"-B";
        argv[2] = (char*)"10001";
        Options opts(argc, argv);

        EXPECT( 10001 == opts.burn() ); 
    },

    CASE( "Options_parameters" )
    {
        int argc=0;
        char **argv;
        Options opts(argc, argv);
        opts.parameters();
    },

    CASE( "Options_parameter_proposal_variance" )
    {
        int argc=3;
        char *argv[3];
        argv[1] = (char*)"-c";
        argv[2] = (char*)"50.50";
        Options opts(argc, argv);
        EXPECT( 50.50 == opts.parameter_proposal_variance() ); 
    },

    CASE( "Options_observation_noise_variance" )
    {
        int argc=3;
        char *argv[3];
        argv[1] = (char*)"-o";
        argv[2] = (char*)"0.4111";
        Options opts(argc, argv);

        EXPECT( 0.4111 == opts.observation_noise_variance() ); 
    },

    CASE( "Options_trajectory_path_delta" )
    {
        int argc=3;
        char *argv[3];
        argv[1] = (char*)"-p";
        argv[2] = (char*)"0.42";
        Options opts(argc, argv);
        EXPECT( 0.42 == opts.trajectory_path_delta() ); 
    },

    CASE( "Options_diffusion_coefficient" )
    {
        int argc=3;
        char *argv[3];
        argv[1] = (char*)"-d";
        argv[2] = (char*)"0.5";
        Options opts(argc, argv);
        EXPECT( 0.5 == opts.diffusion_coefficient() );
    },

};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv /*, std::cout */ );
}
