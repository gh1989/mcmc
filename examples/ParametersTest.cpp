#include "Parameters.h"
#include "LangevinDynamics.h"

using namespace MCMC;

int main()
{
    Parameters< LangevinDynamics::ParameterType, 
                LangevinDynamics::ParameterPointType> params(1);
         
    return 0;
}