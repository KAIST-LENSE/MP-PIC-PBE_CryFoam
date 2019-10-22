/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    mppicPbeFoam

Description
    Transient solver for multiphase particle in cell coupled with populationa balance euqaton (MP-PIC-PBE). The solver is developed by Shin-Hyuk Kim, Jay H. Lee at KAIST and Richard D. Braatz at MIT. The users are refered to a reference for details on this source code; S. H. Kim, R. D. Braatz, and J. H. Lee, "Multiphase particle in cell coupled with population balance equation for multiscale computational fluid dynamcis", Comput. Chem. Eng. 2019.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "rhoReactionThermo.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "crystallizationKinetics.H"

#include "basicKinematicMPPICCloud.H"
#define basicKinematicTypeCloud basicKinematicMPPICCloud

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (!LTS)
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mu = thermo.mu();

        volScalarField x_solute(Y[0]) ;
	volScalarField x_solvent(Y[1]) ;
	volScalarField x_antisolvent(Y[2]) ;

	forAll(mesh.C(), cell)
	{
	    c_sat[cell] = solubility(x_antisolvent[cell], x_solvent[cell], T[cell]);
	    supersat[cell] = (x_solute[cell]/(x_solvent[cell]+x_antisolvent[cell])) - c_sat[cell] ;
	    relsupersat[cell] = (x_solute[cell]/(x_solvent[cell]+x_antisolvent[cell])) / c_sat[cell] ;
            R_growth[cell] = growthRate(supersat[cell], relsupersat[cell]);
	    R_nucleation[cell] = nucRate(supersat[cell], relsupersat[cell]);
	}

        Info<< "Evolving " << kinematicCloud.name() << endl;

        kinematicCloud.evolve();

        // Update continuous phase volume fraction field
        alpha = max(1.0 - kinematicCloud.theta(), alphaMin);
        alpha.correctBoundaryConditions();
        alphaf = fvc::interpolate(alpha);

        fvVectorMatrix cloudSU(kinematicCloud.SU(U));
        volVectorField cloudVolSUSu
        (
            IOobject
            (
                "cloudVolSUSu",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedVector
            (
                "0",
                cloudSU.dimensions()/dimVolume,
                Zero
            ),
            zeroGradientFvPatchVectorField::typeName
        );

        cloudVolSUSu.primitiveFieldRef() = -cloudSU.source()/mesh.V();
        cloudVolSUSu.correctBoundaryConditions();
        cloudSU.source() = Zero;

        R_solute = kinematicCloud.R_total();
        R_solute.correctBoundaryConditions();

        if (pimple.nCorrPIMPLE() <= 1)
        {
            #include "rhoEqn.H"
        }

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "YEqn.H"
            #include "EEqn.H"
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = thermo.rho();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
