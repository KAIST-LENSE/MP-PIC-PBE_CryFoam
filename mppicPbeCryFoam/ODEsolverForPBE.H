#include <math.h>
#include <algorithm>

namespace Foam
{
    double problem(double dr, double R_growth, double R_nucleation, double fj0, double fj1, double fj2, double fj3, double minmodeTheta);
//    double PBERHSRK(double dt, double dr, double R_growth, double R_nucleation, double fj0, double fj1, double fj2, double fj3, double minmodeTheta);
} 

namespace Foam
{
    double problem(double dr, double R_growth, double R_nucleation, double fj0, double fj1, double fj2, double fj3, double minmodeTheta)
    {
        double fr0 = minmod(fj0, fj1, fj2, dr, minmodeTheta);
        double fr1 = minmod(fj1, fj2, fj3, dr, minmodeTheta);
        double result = (-(1/dr) * (R_growth * (fj2 + dr/2 * fr1) - R_growth * (fj1 + dr/2 * fr0) - R_nucleation)); 

        return result; 
    }
/*
    double PBERHSRK(double dt, double dr, double R_growth, double R_nucleation, double fj0, double fj1, double fj2, double fj3, double minmodeTheta)
    {
 	double k1, k2, k3, k4 ;

	k1 = dt*problem(dr, R_growth, R_nucleation, fj0, fj1, fj2, fj3, minmodeTheta);
        k2 = dt*problem(dr, R_growth, R_nucleation, fj0, fj1, fj2 + k1/2, fj3, minmodeTheta);
        k3 = dt*problem(dr, R_growth, R_nucleation, fj0, fj1, fj2 + k2/2, fj3, minmodeTheta);
        k4 = dt*problem(dr, R_growth, R_nucleation, fj0, fj1, fj2 + k3, fj3, minmodeTheta);

        double result = (0.1666667)*(k1+2*k2+2*k3+k4); 

        return result;  
    }
*/
} // End namespace Foam

