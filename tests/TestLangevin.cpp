#include "LangevinDynamics.h"
#include "Options.h"
#include "lest/lest.hpp"

using namespace std;
using namespace MCMC;

const lest::test specification[] =
{

    CASE( "LangevinDynamics_generate_true_trajectory" )
    {   
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.generate_true_trajectory( r );
    },

    CASE( "LangevinDynamics_generate_observations" )
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.generate_observations( r );
    },

    CASE( "LangevinDynamics_propose_trajectory" )
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.propose_trajectory(r);
    },

    CASE( "LangevinDynamics_propose_parameters" )
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.propose_parameters(r);
    },

    CASE( "LangevinDynamics_accept" )
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.accept();
    },

    CASE( "LangevinDynamics_finish" )
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.finish();
    },

    CASE( "LangevinDynamics_store_chain" )
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.store_chain(0);
    },

    CASE( "LangevinDynamics_log_path_ratio" )
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.log_path_ratio();
    },

    CASE( "LangevinDynamics_log_prior_ratio" )
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.log_prior_ratio();
    },
    
    CASE( "LangevinDynamics_log_transition_density_ratio" )
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.log_transition_density_ratio();
    },

    CASE( "LangevinDynamics_propose_parameters" )
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.propose_parameters(r);    
    },

    CASE( "LangevinDynamics_log_transition_density_ratio" )
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.log_transition_density_ratio();    
    },

    CASE( "LangevinDynamics_log_prior_ratio" )
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.log_prior_ratio();    
    },

    CASE( "LangevinDynamics_accept" )
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.accept();
    },

    CASE( "LangevinDynamics_store_chain")
    {
        Options<LangevinDynamics> opts;
        opts.set_parallel_paths(1);
        opts.set_mcmc_trials(1);
        opts.set_path_length(8);
        opts.set_extra_data_ratio(1);
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        LangevinDynamics path_scheme( opts );
        path_scheme.store_chain(0);
    },

};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv , std::cout  );
}
