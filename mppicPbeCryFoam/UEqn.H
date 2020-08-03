fvVectorMatrix UEqn
(
    fvm::ddt(alpha, rho, U) + fvm::div(alpharhophi, U)
  + turbulence->divDevRhoReff(U)
 ==
    cloudSU
  - R_solute*U
  + alpha*rho*g
);

UEqn.relax();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
surfaceScalarField phiForces(rhorAUf*fvc::flux(cloudVolSUSu) + alphaf*rhorAUf*(g & mesh.Sf());


if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            phiForces/rhorAUf - fvc::snGrad(p)*mesh.magSf()
        )
    );


    K = 0.5*magSqr(U);
}
