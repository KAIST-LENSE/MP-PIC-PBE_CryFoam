    rho = thermo.rho();
    rhof = fvc::interpolate(rho);

    const volScalarField psip0(psi*p);
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           fvc::interpolate(rho)*fvc::flux(HbyA)
         + alphaf*rhorAUf*fvc::ddtCorr(rho, U, rhophi)
        )
    );

    if (p.needReference())
    {
        adjustPhi(phiHbyA, U, p);
    }

    phiHbyA += phiForces;

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

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = (phiHbyA + pEqn.flux()/alphaf)/rhof;
            alphaphi = alphaf*phi;
            rhophi = rhof*phi;
            alpharhophi = alphaf*rhof*phi;

	    p.relax();

            U = HbyA + rAU*fvc::reconstruct((phiForces + pEqn.flux()/alphaf)/rhorAUf);
	    U.correctBoundaryConditions();
        }
    }

// Thermodynamic density update
thermo.correctRho(psi*p - psip0);

//#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

K = 0.5*magSqr(U);

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
