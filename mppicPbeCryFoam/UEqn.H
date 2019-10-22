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

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    K = 0.5*magSqr(U);
}
