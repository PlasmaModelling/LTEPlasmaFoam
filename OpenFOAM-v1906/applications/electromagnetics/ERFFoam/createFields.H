	Info<< "Reading transportProperties\n" << endl;

	IOdictionary transportProperties
	(
		IOobject
		(
			"transportProperties",
			runTime.constant(),
			mesh,
			IOobject::MUST_READ_IF_MODIFIED,
			IOobject::NO_WRITE
		)
	);

	dimensionedScalar sigmam
	(
		"sigmam",
		dimensionSet(-1, -3, 3, 0, 0, 2, 0),
		transportProperties
	);

	dimensionedScalar omegam
	(
		"omegam",
		dimensionSet(0, 0, -1, 0, 0, 0, 0),
		transportProperties
	);

	const dimensionedScalar& mu0_ = constant::electromagnetic::mu0;

	Info << "Reading field phi\n" << endl;
	volScalarField phi
	(
		IOobject
		(
			"phi",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	);

// 	label phiRefCell = 0;
// 	scalar phiRefValue = 0.0;
// 	setRefCell(phi, simple.dict(), phiRefCell, phiRefValue);
// 	mesh.setFluxRequired(phi.name());

	Info << "Reading field Jc\n" << endl;
	volVectorField Jc
	(
		IOobject
		(
			"Jc",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		-sigmam*fvc::grad(phi)
	);

	Info<< "Reading field Ar\n" << endl;
	volVectorField Ar
	(
		IOobject
		(
			"Ar",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "Reading field Ai\n" << endl;
	volVectorField Ai
	(
		IOobject
		(
			"Ai",
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
		fvc::curl(Ar - Ai)
	);

	volVectorField E
	(
		IOobject
		(
			"E",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		omegam*(Ar + Ai)
	);

	volVectorField Jind
	(
		IOobject
		(
			"Jind",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		sigmam*E
	);

// 	Info<< nl << endl;
