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

	const dimensionedScalar& mu0_ = constant::electromagnetic::mu0;
	
	Info << "Reading field sigmam\n" << endl;
	volScalarField sigmam
	(
		IOobject
		(
			"sigmam",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		),
		mesh
	);
	
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

// 	label VcRefCell = 0;
// 	scalar VcRefValue = 0.0;
// 	setRefCell(Vc, simple.dict(), VcRefCell, VcRefValue);
// 	mesh.setFluxRequired(Vc.name());

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
		sigmam*E
	);

	Info<< nl << endl;
