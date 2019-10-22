#include <math.h>
#include <algorithm>

namespace Foam
{
    double growthRate(double supersat, double relSupersat);
    double nucRate(double supersat, double relSupersat);
    double solubility(double antiSolvent, double solvent, double T);
} 

namespace Foam
{
    double growthRate(double supersat, double relSupersat)
    {
        double G;

        if(supersat > 0.0)
        {
            G = 8.333e-30*pow((2.4623e3*log(relSupersat)),6.7);
        }
        else
        {
            G = 0.0;
        }

        return G;
    }

    double nucRate(double supersat, double relSupersat)
    {
        double B_homo, B_hetero, B;

        if(supersat > 0.0)
        {
            B_homo = 6.9656e14 * exp(-15.7647/(log(relSupersat) * log(relSupersat)));
            B_hetero = 2.19196e8 * exp(-0.994387/(log(relSupersat) * log(relSupersat)));
            B = B_homo + B_hetero;
        }
        else
        {
            B_homo = 0.0;
            B_hetero = 0.0;
            B = 0.0;
        }

        return B;
    }

    double solubility(double antiSolvent, double solvent, double T)
    {
        double csat, Was, theta;

        if((antiSolvent+solvent)>0.0)
        {
            Was = 100.0*antiSolvent/(antiSolvent+solvent);
        }
        else
        {
            Was = 0.0;
        }

        Was = std::min(std::max(Was, 0.0),100.0);
        theta = T/296.0;
        if(Was <= 45.67)
        {
            csat = 0.001*exp(15.45763*(1.0-1.0/theta))*(-2.7455e-4*pow(Was,3.0)+3.3716e-2*pow(Was,2)-1.6704*Was+33.089);
        }
        else 
        {
            csat = 0.001*exp(15.45763*(1.0-1.0/theta))*(-1.7884e-2*Was+1.7888);
        }

        return csat;
    }
} // End namespace Foam

