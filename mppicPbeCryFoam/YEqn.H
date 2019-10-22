tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        alpharhophi,
        mesh.divScheme("div(alpharhophi,Yi_h)")
    )
);

{
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];
            volScalarField s = R_solute ;
            if(i != 0)
            {
                s = 0 * R_solute ;
            }

            fvScalarMatrix YEqn
            (
                fvm::ddt(alpha, rho, Yi)
              + mvConvection->fvmDiv(alpharhophi, Yi)
             ==
                fvm::laplacian(turbulence->muEff(), Yi)
              - s
            );

            YEqn.relax();

            YEqn.solve(mesh.solver("Yi"));

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
