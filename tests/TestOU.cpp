#include "OUDynamics.h"
#include "OUDynamics.h"
#include "Options.h"
#include "lest/lest.hpp"

using namespace std;
using namespace MCMC;

const lest::test specification[] =
{

    CASE( "OUDynamics_generate_true_trajectory" )
    {   
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics path_scheme( opts );
        path_scheme.generate_true_trajectory( r );
    },

    CASE( "OUDynamics_generate_observations" )
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics path_scheme( opts );
        path_scheme.generate_observations( r );
    },

    CASE( "OUDynamics_propose_trajectory" )
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics path_scheme( opts );
        path_scheme.propose_trajectory(r);
    },

    CASE( "OUDynamics_propose_parameters" )
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics path_scheme( opts );
        path_scheme.propose_parameters(r);
    },

    CASE( "OUDynamics_accept" )
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics path_scheme( opts );
        path_scheme.accept();
    },

    CASE( "OUDynamics_finish" )
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics path_scheme( opts );
        path_scheme.finish();
    },

    CASE( "OUDynamics_store_chain" )
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics path_scheme( opts );
        path_scheme.store_chain(0);
    },

    CASE( "OUDynamics_log_path_ratio" )
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics path_scheme( opts );
        path_scheme.log_path_ratio();
    },

    CASE( "OUDynamics_log_prior_ratio" )
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics path_scheme( opts );
        path_scheme.log_prior_ratio();
    },
    
    CASE( "OUDynamics_log_transition_density_ratio" )
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics path_scheme( opts );
        path_scheme.log_transition_density_ratio();
    },

    CASE( "OUDynamics_propose_parameters" )
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics params( opts );
        params.propose_parameters(r);    
    },

    CASE( "OUDynamics_log_transition_density_ratio" )
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics params( opts );
        params.log_transition_density_ratio();    
    },

    CASE( "OUDynamics_log_prior_ratio" )
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics params( opts );
        params.log_prior_ratio();    
    },

    CASE( "OUDynamics_accept" )
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics params( opts );
        params.accept();
    },

    CASE( "OUDynamics_store_chain")
    {
        Options<OUDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        OUDynamics params( opts );
        params.store_chain(0);
    },

};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv , std::cout  );
}
