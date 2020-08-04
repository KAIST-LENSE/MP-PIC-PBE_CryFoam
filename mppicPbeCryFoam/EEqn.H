{
    dHmix = (-5.90806E+05*pow(x_antisolvent,6.0)+1.45022E+06*pow(x_antisolvent,5.0)-9.73106E+05*pow(x_antisolvent,4.0)+3.90128E+04*pow(x_antisolvent,3.0)+2.27166E+05*pow(x_antisolvent,2.0)-1.52557E+05*x_antisolvent)*dHmixDimensions;
    
    forAll(mesh.C(), cell)
    {
        dHmixCell[cell] = dHmix[cell];
    }
 
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(alpha, rho, he) + fvm::div(alpharhophi, he)
      + alpha*rho*(dHmixCell-dHmixCell.oldTime().oldTime())/runTime.deltaT()
      + alpha*dHmixCell*fvc::ddt(rho)
      + rho*dHmixCell*fvc::ddt(alpha)
      + fvc::div(alpharhophi, dHmix)
      + fvc::ddt(alpha, rho, K) + fvc::div(alpharhophi, K)
      - fvm::laplacian(alpha*turbulence->alphaEff(), he)
      + (
            he.name() == "e"
          ? fvc::div(fvc::absolute(alphaphi, U), p, "div(phiv,p)")
          : -alpha*dpdt
        )
    ==
        alpha*rho*(U&g)
      + R_solute*heatOfCrystal
    );
    EEqn.relax();

    EEqn.solve();

    thermo.correct();

    Info<< "T min/max   " << min(T).value() << ", "
        << max(T).value() << endl;
}
