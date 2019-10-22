if (!pimple.SIMPLErho())
{
    rho = thermo.rho();
}
    rhof = fvc::interpolate(rho);
    const volScalarField psip0(psi*p);
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           fvc::interpolate(rho)*fvc::flux(HbyA)
         + alphaf*rhorAUf*fvc::ddtCorr(rho, U, rhophi)
        )
    );

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAUf);

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(alpha, rho)
      + alpha*psi*correction(fvm::ddt(p))
      + fvc::div(alphaf*phiHbyA)
      ==
      - R_solute
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(alphaf*rhorAUf, p));

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = (phiHbyA + pEqn.flux()/alphaf)/rhof;
            alphaphi = alphaf*phi;
            rhophi = rhof*phi;
            alpharhophi = alphaf*rhof*phi;
        }
    }

// Thermodynamic density update
thermo.correctRho(psi*p - psip0);

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

p.relax();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
K = 0.5*magSqr(U);


p.correctBoundaryConditions();

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
/*
rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
Info<< "rho min/max = " << min(rho).value() << ", " << max(rho).value() << endl;
*/
