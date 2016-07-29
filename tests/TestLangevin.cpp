#include "LangevinPathScheme.h"
#include "LangevinParams.h"
#include "Options.h"
#include "lest/lest.hpp"

using namespace std;
using namespace MCMC;

const lest::test specification[] =
{

    CASE( "LangevinPathScheme_generate_true_trajectory" )
    {   
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinPathScheme path_scheme( opts );
        path_scheme.generate_true_trajectory( r );
    },

    CASE( "LangevinPathScheme_generate_observations" )
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinPathScheme path_scheme( opts );
        path_scheme.generate_observations( r );
    },

    CASE( "LangevinPathScheme_propose_trajectory" )
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinPathScheme path_scheme( opts );
        path_scheme.propose_trajectory(r);
    },

    CASE( "LangevinPathScheme_propose_parameters" )
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinPathScheme path_scheme( opts );
        path_scheme.propose_parameters(r);
    },

    CASE( "LangevinPathScheme_accept" )
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinPathScheme path_scheme( opts );
        path_scheme.accept();
    },

    CASE( "LangevinPathScheme_finish" )
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinPathScheme path_scheme( opts );
        path_scheme.finish();
    },

    CASE( "LangevinPathScheme_store_chain" )
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinPathScheme path_scheme( opts );
        path_scheme.store_chain(0);
    },

    CASE( "LangevinPathScheme_log_path_ratio" )
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinPathScheme path_scheme( opts );
        path_scheme.log_path_ratio();
    },

    CASE( "LangevinPathScheme_log_prior_ratio" )
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinPathScheme path_scheme( opts );
        path_scheme.log_prior_ratio();
    },
    
    CASE( "LangevinPathScheme_log_transition_density_ratio" )
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinPathScheme path_scheme( opts );
        path_scheme.log_transition_density_ratio();
    },

    CASE( "LangevinParams_propose_parameters" )
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinParams params( opts );
        params.propose_parameters(r);    
    },

    CASE( "LangevinParams_log_transition_density_ratio" )
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinParams params( opts );
        params.log_transition_density_ratio();    
    },

    CASE( "LangevinParams_log_prior_ratio" )
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinParams params( opts );
        params.log_prior_ratio();    
    },

    CASE( "LangevinParams_accept" )
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinParams params( opts );
        params.accept();
    },

    CASE( "LangevinParams_store_chain")
    {
        Options opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinParams params( opts );
        params.store_chain(0);
    },

};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv , std::cout  );
}
