#include <math.h>
#include <algorithm>

namespace Foam
{
    double minmod(double fL, double f, double fR, double dr, double theta);
} 

namespace Foam
{
    double minmod(double fL, double f, double fR, double dr, double theta)
    {
        double grad1, grad2, grad3, grad;

        grad1 = theta*(f-fL)/dr;
        grad2 = 0.5*(fR-fL)/dr;
        grad3 = theta*(fR-f)/dr;

        if ((grad1>0.0) && (grad2>0.0) && (grad3>0.0)) 
        {
            grad = std::min(grad1,grad2);
            grad = std::min(grad,grad3);
        }
        else if ((grad1<0.0) && (grad2<0.0) && (grad3<0.0)) 
        {
            grad = std::max(grad1,grad2);
            grad = std::max(grad, grad3);
        }
        else
        {
            grad = 0.0;
        }

        return grad;  
    }
} // End namespace Foam

