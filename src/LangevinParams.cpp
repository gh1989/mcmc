#include "LangevinParams.h"
#include <iostream>
#include <fstream>

using namespace MCMC;

void LangevinParamImp::propose_parameters(gsl_rng *r)
{
    double variance = opts.parameter_proposal_variance();
 
    for( size_t i=0; i<defining_modes; ++i)
    {
        c_star(i) = ComplexType(gsl_ran_gaussian(r, variance ), gsl_ran_gaussian(r, variance)) + c(i);
    }
    update_potential_star();
}

void LangevinParamImp::update_potential_star()
{
    int idx;

    int ii;
    int jj; 
 
    for( size_t i=1; i<cutoff+1; ++i )
        _potential_star.set_mode( i, 0, c_star(i-1) );

    for( size_t j=1; j<cutoff+1; ++j )
        for( signed int i=-cutoff; i<cutoff+1; ++i )
        { 
            idx = cutoff + (j-1)*(2*cutoff+1) + (i+cutoff);
            _potential_star.set_mode( i, j, c_star( idx ) ); 
        }

}

void LangevinParamImp::update_potential()
{
    int idx;

    int ii;
    int jj; 
   
    for( size_t i=1; i<cutoff+1; ++i )
        _potential.set_mode( i, 0, c(i-1) );

    for( int i=-cutoff; i<cutoff+1; ++i )
        for( size_t j=1; j<cutoff+1; ++j ){ 
            idx = cutoff + (j-1)*(2*cutoff+1) + (i+cutoff);
            _potential.set_mode( i, j, c( idx ) );
        }

}

double LangevinParamImp::log_transition_density_ratio()
{
    // Symmetric
    return 0;
}

double LangevinParamImp::log_prior_ratio()
{
    double log_total = 0;
    double var = opts.parameter_proposal_variance();
    double normal_constant = 0.5 / var;

    for( size_t i=0; i<defining_modes; ++i)
        log_total -= normal_constant*( pow( std::abs(c_star(i)), 2) -  pow( std::abs(c(i)), 2) );

    return log_total;
}

void LangevinParamImp::accept()
{
    for( size_t i=0; i<defining_modes; ++i ) 
        c(i) = c_star(i);
    update_potential();
    _potential.print_modes();
}

void LangevinParamImp::store_chain(int n)
{
    for( size_t i=0; i<defining_modes; ++i)
    {
        constants_chain(n,i) = c(i);
    }
}

void LangevinParamImp::output_file_timeseries()
{
    size_t N  = opts.mcmc_trials();

    std::ofstream mcmc_file;
    mcmc_file.open ("output/langevin_mcmc_timeseries.txt");

    for( size_t n=0; n<N; ++n )
    {
        for( size_t m=0; m<defining_modes; ++m)
            mcmc_file << std::real(constants_chain(n,m)) << "\t" << std::imag(constants_chain(n,m)) << "\t";

        mcmc_file << std::endl;
    }
    mcmc_file.close();
} 
