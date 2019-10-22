#include "createRDeltaT.H"

#include "readGravitationalAcceleration.H"

Info<< "Creating themophysical properties model\n" << endl;

autoPtr<rhoReactionThermo> Thermo
(
    rhoReactionThermo::New(mesh)
);

rhoReactionThermo& thermo = Thermo();

thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

volScalarField& p = thermo.p();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField c_sat
(
    IOobject
    (
        "c_sat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("c_sat", dimensionSet (0,0,0,0,0,0,0), 0.0)
);

volScalarField supersat
(
    IOobject
    (
        "supersat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("supersat", dimensionSet (0,0,0,0,0,0,0), 0.0)
);

volScalarField relsupersat
(
    IOobject
    (
        "relsupersat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("relsupersat", dimensionSet (0,0,0,0,0,0,0), 0.0)
);

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.mu()
);

volScalarField R_growth
(
    IOobject
    (
        "R_growth",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_growth", dimensionSet (0,1,-1,0,0,0,0), 0.0)
);

volScalarField R_nucleation
(
    IOobject
    (
        "R_nucleation",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_nucleation", dimensionSet (0,-3,-1,0,0,0,0), 0.0)
);

Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

dimensionedScalar rhoMax
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMax",
        pimple.dict(),
        dimDensity,
        GREAT
    )
);

dimensionedScalar rhoMin
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMin",
        pimple.dict(),
        dimDensity,
        0
    )
);

Info<< "Creating multi-variate interpolation scheme\n" << endl;
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
);

volScalarField dHmix
(
  IOobject
  (
    "dHmix",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("dHmix", dimensionSet (0, 2,-2,0,0,0,0), 0.0)
);

volScalarField dHmixCell
(
  IOobject
  (
    "dHmixCell",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("dHmixCell", dimensionSet (0, 2,-2,0,0,0,0), 0.0)
);

volScalarField dHmixDimensions
(
  IOobject
  (
    "dHmixDimensions",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("dHmixDimensions", dimensionSet (0, 2,-2,0,0,0,0), 1.0)
);

/*
volScalarField H_realAnti
(
  IOobject
  (
    "H_realAnti",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("H_realAnti", dimensionSet (0,2,-2,0,0,0,0), 0.0)
);

volScalarField H_realSolvent
(
  IOobject
  (
    "H_realSolvent",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("H_realSolvent", dimensionSet (0,2,-2,0,0,0,0), 0.0)
);
*/
Info << "Creating field alphac\n" << endl;
// alphac must be constructed before the cloud
// so that the drag-models can find it
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);

volScalarField R_solute
(
    IOobject
    (
        "R_solute",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_solute", dimensionSet (1,-3,-1,0,0,0,0), 0.0)
);

volScalarField heatOfCrystal
(
    IOobject
    (
        "heatOfCrystal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("heatOfCrystal", dimensionSet (0, 2,-2,0,0,0,0), 38042500/404.256)
);

#include "createClouds.H"

// Particle fraction upper limit
scalar alphaMin
(
    1.0
  - readScalar
    (
        kinematicCloud.particleProperties().subDict("constantProperties")
       .lookup("alphaMax")
    )
);

// Update alphac from the particle locations
alpha = max(1.0 - kinematicCloud.theta(), alphaMin);
alpha.correctBoundaryConditions();
R_solute = kinematicCloud.R_total() ;
R_solute.correctBoundaryConditions();

surfaceScalarField alphaf("alphaf", fvc::interpolate(alpha));
surfaceScalarField rhof("rhof", fvc::interpolate(rho));
surfaceScalarField alphaphi("alphaphi", alphaf*phi);
surfaceScalarField rhophi("rhophi", rhof*phi);
surfaceScalarField alpharhophi("alpharhophi", alphaf*rhof*phi);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        alpharhophi,
        thermo
    )
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));



