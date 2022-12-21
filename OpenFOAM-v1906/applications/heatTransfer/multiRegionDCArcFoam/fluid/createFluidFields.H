// Initialise fluid field pointer lists
PtrList<rhoThermo> thermoFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volScalarField> muFluid(fluidRegions.size());
PtrList<volScalarField> alphaFluid(fluidRegions.size());
PtrList<volScalarField> CpFluid(fluidRegions.size());
PtrList<volScalarField> sigmaEFluid(fluidRegions.size());
PtrList<volScalarField> VcFluid(fluidRegions.size());
PtrList<volVectorField> AFluid(fluidRegions.size());
PtrList<volVectorField> BFluid(fluidRegions.size());
PtrList<volVectorField> EcFluid(fluidRegions.size());
PtrList<volVectorField> JFluid(fluidRegions.size());
PtrList<volVectorField> LoFluid(fluidRegions.size());
PtrList<volScalarField> JoFluid(fluidRegions.size());
PtrList<volScalarField> EJFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<compressible::turbulenceModel> turbulence(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<radiation::radiationModel> radiation(fluidRegions.size());

List<scalar> initialMassFluid(fluidRegions.size());
List<label> pRefCellFluid(fluidRegions.size(), -1);
List<scalar> pRefValueFluid(fluidRegions.size(), Zero);
List<bool> frozenFlowFluid(fluidRegions.size(), false);

PtrList<dimensionedScalar> rhoMax(fluidRegions.size());
PtrList<dimensionedScalar> rhoMin(fluidRegions.size());

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoFluid\n" << endl;

    thermoFluid.set
    (
        i,
        rhoThermo::New(fluidRegions[i]).ptr()
    );

	const dimensionedScalar& mu0_ = constant::electromagnetic::mu0;
	const dimensionedScalar& e_ = constant::electromagnetic::e;
	const dimensionedScalar& kb_ = constant::physicoChemical::k;

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].rho()
        )
    );
	
	muFluid.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"mu",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			thermoFluid[i].mu()
		)
	);
	
	alphaFluid.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"alpha",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			thermoFluid[i].alpha()
		)
	);

	CpFluid.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"Cp",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			thermoFluid[i].Cp()
		)
	);

	sigmaEFluid.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"sigmaE",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		)
	);
// 	sigmaE = thermoFluid[i].sigmaE();
// 	sigmaE.correctBoundaryConditions();
	
	Info << "Reading fields Vc\n" << endl;
	VcFluid.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"Vc",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		)
	);
	
	Info<< "Reading fields A\n" << endl;
	AFluid.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"A",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		)
	);
	
	BFluid.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"B",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			fvc::curl(AFluid[i])
		)
	);
	
	EcFluid.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"Ec",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			-fvc::grad(VcFluid[i])
		)
	);
	
	JFluid.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"J",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			sigmaEFluid[i]*EcFluid[i]
		)
	);
	
	LoFluid.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"Lo",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			JFluid[i]^BFluid[i]
		)
	);
	
	JoFluid.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"Jo",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			magSqr(JFluid[i])/sigmaEFluid[i]
		)
	);
	
	EJFluid.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"EJ",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			fvc::div(5*kb_*JFluid[i]*thermoFluid[i].he()/(2*e_*CpFluid[i]))
		)
	);

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoFluid[i]*UFluid[i])
                & fluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, Zero) // uses name
        )
    );

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding to turbulence\n" << endl;
    turbulence.set
    (
        i,
        compressible::turbulenceModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // Force p_rgh to be consistent with p
    p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    radiation.set
    (
        i,
        radiation::radiationModel::New(thermoFluid[i].T())
    );

    initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

    const dictionary& simpleDict =
        fluidRegions[i].solutionDict().subDict("SIMPLE");

    setRefCell
    (
        thermoFluid[i].p(),
        p_rghFluid[i],
        simpleDict,
        pRefCellFluid[i],
        pRefValueFluid[i]
    );

    simpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

    rhoMax.set
    (
        i,
        new dimensionedScalar("rhoMax", dimDensity, GREAT, simpleDict)
    );

    rhoMin.set
    (
        i,
        new dimensionedScalar("rhoMin", dimDensity, Zero, simpleDict)
    );

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );

    turbulence[i].validate();
}
