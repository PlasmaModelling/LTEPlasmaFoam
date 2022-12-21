Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));
rhoThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

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

volScalarField alpha
(
	IOobject
	(
		"alpha",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	thermo.alpha()
);

volScalarField Cp
(
	IOobject
	(
		"Cp",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	thermo.Cp()
);

volScalarField& he = thermo.he();

volScalarField Energy
(
	IOobject
	(
		"Energy",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	he
);

volScalarField& p = thermo.p();

const dimensionedScalar& mu0_ = constant::electromagnetic::mu0;
const dimensionedScalar& e_ = constant::electromagnetic::e;
const dimensionedScalar& kb_ = constant::physicoChemical::k;
// const dimensionedScalar ome ("ome", dimensionSet(0, 0, -1, 0, 0, 0, 0), 1);

volScalarField sigmaE
(
	IOobject
	(
		"sigmaE",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);
sigmaE = thermo.sigmaE();
sigmaE.correctBoundaryConditions();

Info << "Reading field Vc\n" << endl;
volScalarField Vc
(
	IOobject
	(
		"Vc",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

Info<< "Reading field A\n" << endl;
volVectorField A
(
	IOobject
	(
		"A",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

volVectorField B
(
	IOobject
	(
		"B",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fvc::curl(A)
);

volVectorField Ec
(
	IOobject
	(
		"Ec",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	-fvc::grad(Vc)
);

volVectorField J
(
	IOobject
	(
		"J",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	sigmaE*Ec
);

volVectorField Lo
(
	IOobject
	(
		"Lo",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	J^B
);

// volScalarField dumpF
// (
// 	IOobject
// 	(
// 		"dumpF",
// 		runTime.timeName(),
// 		mesh,
// 		IOobject::NO_READ,
// 		IOobject::NO_WRITE
// 	),
// 	he/he
// );

volScalarField Jo
(
	IOobject
	(
		"Jo",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	magSqr(J)/sigmaE
);

volScalarField EJ //volVectorField EJ
(
	IOobject
	(
		"EJ",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
// 	5*kb_*J*he/(2*e_*Cp)
	fvc::div(5*kb_*J*he/(2*e_*Cp))
);

Info<< "Reading field U\n" << endl;
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

#include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
p_rgh = p - rho*gh;

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    simple.dict(),
    pRefCell,
    pRefValue
);

mesh.setFluxRequired(p_rgh.name());

dimensionedScalar initialMass = fvc::domainIntegrate(rho);
dimensionedScalar totalVolume = sum(mesh.V());

#include "createMRF.H"
#include "createRadiationModel.H"

const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, simple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, simple.dict());

#include "createFvOptions.H"
