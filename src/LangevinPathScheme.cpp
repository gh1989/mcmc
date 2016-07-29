#include "LangevinPathScheme.h"
#include <iostream>

using namespace MCMC;

double LangevinPathSchemeImp::log_path_ratio()
{
    FourierSeries v = params().potential();
    FourierSeries v_star = params().potential_star();

    // Diffusion sigma
    double d_sigma = opts.diffusion_coefficient();  
    double dt = opts.trajectory_path_delta(); 
    double diff_const = (0.5 / (d_sigma*d_sigma*dt));

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
            log_total -= obs_const * pow( proposed(k, l, 0, 0 ) - observed(k, l, 0), 2);
            log_total -= obs_const * pow( proposed(k, l, 0, 1 ) - observed(k, l, 1), 2);

            log_total += obs_const * pow( current(k, l, 0, 0 ) - observed(k, l, 0), 2);
            log_total += obs_const * pow( current(k, l, 0, 1 ) - observed(k, l, 1), 2);

            for( size_t m=1; m<M; ++m )
            {
                xt << proposed(k, l, m, 0 ), proposed(k, l, m, 1 );
                xtminus1 << proposed(k, l, m-1, 0 ), proposed(k, l, m-1, 1 );
                gradvtminus1 = v_star.grad( xtminus1 ); 

                log_total -= diff_const * pow( xt(0)-xtminus1(0)+gradvtminus1(0)*dt, 2);
                log_total -= diff_const * pow( xt(0)-xtminus1(1)+gradvtminus1(1)*dt, 2);

                xt << current(k, l, m, 0 ), current(k, l, m, 1 );
                xtminus1 << current(k, l, m-1, 0 ), current(k, l, m-1, 1 );
                gradvtminus1 = v.grad( xtminus1 ); 

                log_total += diff_const * pow( xt(0)-xtminus1(0)+gradvtminus1(0)*dt, 2);
                log_total += diff_const * pow( xt(0)-xtminus1(1)+gradvtminus1(1)*dt, 2);

            }
        }
    return log_total;
}

void LangevinPathSchemeImp::generate_true_trajectory(gsl_rng *r)
{
    generate_random_starts( r, real );
    langevin_trajectory(r, true_potential, real );
}

void LangevinPathSchemeImp::generate_observations(gsl_rng *r)
{
    double variance = opts.observation_noise_variance();

    double random_noise_x;
    double random_noise_y;

    size_t K = opts.parallel_paths();
    size_t L = opts.path_length();

    for( size_t k=0; k<K; ++k )
    {
        for( size_t l=0; l<L; ++l )
        {   
            random_noise_y = gsl_ran_gaussian(r, variance );
            random_noise_x = gsl_ran_gaussian(r, variance );                

            observed(k, l, 0 ) = real(k, l, 0, 0 ) + random_noise_x;
            observed(k, l, 1 ) = real(k, l, 0, 1 ) + random_noise_y;
        }   
    }
}

void LangevinPathSchemeImp::propose_trajectory(gsl_rng *r)
{
    FourierSeries v_star = params().potential_star();
    setup_observed_starts( r, observed, proposed );
    langevin_trajectory(r, v_star, proposed );
}

///////////////////////////////////////////////////////////////////////////////

void LangevinPathSchemeImp::generate_random_starts( gsl_rng *r, Tensor<double, 4> &out )
{
    size_t K = opts.parallel_paths();

    for( size_t k=0; k<K; ++k )
    {
        // Initialisation: random point inside [-1,1]X[-1,1].
        out(k, 0, 0, 0 ) = 2.0*gsl_rng_uniform(r)-1;
        out(k, 0, 0, 1 ) = 2.0*gsl_rng_uniform(r)-1;
    }
}

void LangevinPathSchemeImp::setup_observed_starts( gsl_rng *r, const Tensor<double, 3> &y, Tensor<double, 4> &out )
{
    double variance = opts.observation_noise_variance(); 
    size_t K = opts.parallel_paths();

    for( size_t k=0; k<K; ++k )
    {
        out(k, 0, 0, 0 ) = y(k, 0, 0); // + gsl_ran_gaussian(r, variance);
        out(k, 0, 0, 1 ) = y(k, 0, 1); // + gsl_ran_gaussian(r, variance);
    }
}

void LangevinPathSchemeImp::langevin_trajectory( gsl_rng *r, FourierSeries &V, Tensor<double, 4> &out )
{
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

