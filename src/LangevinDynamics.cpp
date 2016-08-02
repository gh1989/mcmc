#include "LangevinDynamics.h"
#include <iostream>
#include <fstream>

ComplexType LangevinDynamics::sample_transition_density(gsl_rng *r, ComplexType c )
{
    double variance = opts.parameter_proposal_variance();
    return ComplexType(gsl_ran_gaussian(r, variance ), gsl_ran_gaussian(r, variance)) + c;
}

double LangevinDynamics::log_transition_density_ratio(ParameterType &c_star, ParameterType &c)
{
    // Symmetric
    return 0;
}

double LangevinDynamics::log_prior_ratio(ParameterType &c_star, ParameterType &c)
{
    double log_total = 0;
    double var = opts.parameter_proposal_variance();
    double normal_constant = 0.5 / var;

    int c_dim = parameter_dimension();

    for( size_t i=0; i<c_dim; ++i)
        log_total -= normal_constant*( pow( std::abs(c_star(i)), 2) - pow( std::abs(c(i)), 2) );

    return log_total;
}

double LangevinDynamics::log_path_likelihood(   PathType &x, 
                                                ParameterType &c, 
                                                CoarsePathType &y)
{
    FourierSeries V(_cutoff);
    V.set_modes(c);

    // Diffusion sigma
    double sigma = opts.diffusion_coefficient();  
    double dt = opts.trajectory_path_delta(); 
    double diff_const = (0.5 / (sigma*sigma*dt));

    // Observation variance
    double o_variance = opts.observation_noise_variance();
    double obs_const = (0.5 / o_variance);

    size_t K = opts.parallel_paths();
    size_t M = opts.extra_data_ratio();
    size_t L = opts.path_length(); 
    
    double log_total = 0;

    Vector2d xt;
    Vector2d xtminus1;
    Vector2d gradvtminus1;

    for( size_t k=0; k<K; ++k )
        for( size_t l=0; l<L; ++l )
        {
            log_total -= obs_const * pow( x(k, l, 0, 0 ) - y(k, l, 0), 2);
            log_total -= obs_const * pow( x(k, l, 0, 1 ) - y(k, l, 1), 2);
            
            for( size_t m=1; m<M; ++m )
            {
                xt << x(k, l, m, 0 ), x(k, l, m, 1 );
                xtminus1 << x(k, l, m-1, 0 ), x(k, l, m-1, 1 );
                gradvtminus1 = V.grad( xtminus1 ); 

                log_total -= diff_const * pow( xt(0)-xtminus1(0)+gradvtminus1(0)*dt, 2);
                log_total -= diff_const * pow( xt(1)-xtminus1(1)+gradvtminus1(1)*dt, 2);

            }
        }
    return log_total;
}

void LangevinDynamics::generate_random_starts( gsl_rng *r, PathType &out )
{
    size_t K = opts.parallel_paths();

    for( size_t k=0; k<K; ++k )
    {
        // Initialisation: random point inside [-1,1]X[-1,1].
        out(k, 0, 0, 0 ) = 2.0*gsl_rng_uniform(r)-1;
        out(k, 0, 0, 1 ) = 2.0*gsl_rng_uniform(r)-1;
    }
}

void LangevinDynamics::setup_observed_starts( gsl_rng *r, CoarsePathType &y, PathType &out )
{
    double variance = opts.observation_noise_variance(); 
    size_t K = opts.parallel_paths();

    for( size_t k=0; k<K; ++k )
    {
        out(k, 0, 0, 0 ) = y(k, 0, 0); // + gsl_ran_gaussian(r, variance);
        out(k, 0, 0, 1 ) = y(k, 0, 1); // + gsl_ran_gaussian(r, variance);
    }
}

void LangevinDynamics::trajectory( gsl_rng *r, ParameterType &c, PathType &out )
{
    
    FourierSeries V(_cutoff);
    V.set_modes( c );    
        
    double dt = opts.trajectory_path_delta();
    double sigma = opts.diffusion_coefficient();    

    double random_noise_x;
    double random_noise_y;

    size_t K = opts.parallel_paths();
    size_t M = opts.extra_data_ratio();
    size_t L = opts.path_length();   

    Vector2d xtminus1;
    Vector2d gradtminus1;

    for( size_t k=0; k<K; ++k )
        for( size_t l=0; l<L-1; ++l )
        {  
            for( size_t m=1; m<M; ++m )
            {
                random_noise_x = gsl_ran_gaussian(r, sigma*sigma*dt);
                random_noise_y = gsl_ran_gaussian(r, sigma*sigma*dt);                

                xtminus1 << out(k, l, m-1, 0 ), out(k, l, m-1, 1 );
                gradtminus1 = V.grad( xtminus1 );
        
                out(k, l, m, 0 ) = xtminus1(0)-gradtminus1(0)*dt + random_noise_x;
                out(k, l, m, 1 ) = xtminus1(1)-gradtminus1(1)*dt + random_noise_y;
            }   

            random_noise_x = gsl_ran_gaussian(r, sigma*sigma*dt);
            random_noise_y = gsl_ran_gaussian(r, sigma*sigma*dt);

            xtminus1 << out(k, l, M-1, 0 ), out(k, l, M-1, 1 );
            gradtminus1 = V.grad( xtminus1 );
                
            out(k, l+1, 0, 0) = xtminus1(0)-gradtminus1(0)*dt + random_noise_x;
            out(k, l+1, 0, 1) = xtminus1(1)-gradtminus1(1)*dt + random_noise_y;

        }   
}

void LangevinDynamics::output_file_timeseries(ParameterChainType &ccc)
{
    std::cout<<"LangevinDynamics::output_file_timeseries"<<std::endl;
    size_t N  = opts.mcmc_trials();

    std::ofstream mcmc_file;
    mcmc_file.open ("output/langevin_mcmc_timeseries.txt");

    int defining_modes = parameter_dimension(opts);

    for( size_t n=0; n<N; ++n )
    {
        for( size_t m=0; m<defining_modes; ++m)
            mcmc_file << std::real(ccc(n,m)) << "\t" << std::imag(ccc(n,m)) << "\t";

        mcmc_file << std::endl;
    }
    mcmc_file.close();
} 
