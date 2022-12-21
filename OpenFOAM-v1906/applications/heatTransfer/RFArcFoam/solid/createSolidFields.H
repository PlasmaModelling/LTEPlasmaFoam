    // Initialise solid field pointer lists
    PtrList<coordinateSystem> coordinates(solidRegions.size());
    PtrList<solidThermo> thermos(solidRegions.size());
    PtrList<radiation::radiationModel> radiations(solidRegions.size());
    PtrList<fv::options> solidHeatSources(solidRegions.size());
    PtrList<volScalarField> betavSolid(solidRegions.size());
	PtrList<volScalarField> sigmaSSolid(solidRegions.size());
	PtrList<volScalarField> VcSolid(solidRegions.size());
	PtrList<volVectorField> ArSolid(solidRegions.size());
	PtrList<volVectorField> AiSolid(solidRegions.size());
// 	PtrList<volVectorField> BSolid(solidRegions.size());
// 	PtrList<volVectorField> EcSolid(solidRegions.size());
// 	PtrList<volVectorField> JSolid(solidRegions.size());
	PtrList<volVectorField> JcoilSolid(solidRegions.size());
    PtrList<volSymmTensorField> aniAlphas(solidRegions.size());

    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermos\n" << endl;
        thermos.set(i, solidThermo::New(solidRegions[i]));

        Info<< "    Adding to radiations\n" << endl;
        radiations.set(i, radiation::radiationModel::New(thermos[i].T()));

        Info<< "    Adding fvOptions\n" << endl;
        solidHeatSources.set
        (
            i,
            new fv::options(solidRegions[i])
        );

	IOdictionary thermophysicalDict
	(
		IOobject
		(
			"thermophysicalProperties",    // dictionary name
			runTime.constant(),     // dict is found in "constant"
			solidRegions[i],                   // registry for the dict
			IOobject::MUST_READ,    // must exist, otherwise failure
			IOobject::NO_WRITE      // dict is only read by the solver
		)
	);
	const dictionary& omegaDict = thermophysicalDict.subDict("frequency");

	const dimensionedScalar omegaE
	(
		"omegaE",
		dimensionSet(0, 0, -1, 0, 0, 0, 0),
		readScalar(omegaDict.lookup("omegaE"))
	);

        if (!thermos[i].isotropic())
        {
            Info<< "    Adding coordinateSystems\n" << endl;
            coordinates.set
            (
                i,
                coordinateSystem::New
                (
                    solidRegions[i],
                    thermos[i],
                    coordinateSystem::typeName_()
                )
            );

            tmp<volVectorField> tkappaByCp =
                thermos[i].Kappa()/thermos[i].Cp();

            aniAlphas.set
            (
                i,
                new volSymmTensorField
                (
                    IOobject
                    (
                        "Anialpha",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedSymmTensor(tkappaByCp().dimensions(), Zero),
                    zeroGradientFvPatchSymmTensorField::typeName
                )
            );

            aniAlphas[i].primitiveFieldRef() =
                coordinates[i].transformPrincipal
                (
                    solidRegions[i].cellCentres(),
                    tkappaByCp()
                );
            aniAlphas[i].correctBoundaryConditions();

        }

        IOobject betavSolidIO
        (
            "betavSolid",
            runTime.timeName(),
            solidRegions[i],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (betavSolidIO.typeHeaderOk<volScalarField>(true))
        {
            betavSolid.set
            (
                i,
                new volScalarField(betavSolidIO, solidRegions[i])
            );
        }
        else
        {
            betavSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "betavSolid",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("1", dimless, scalar(1))
                )
            );
        }

			sigmaSSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "sigmaE",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
//                     solidRegions[i]
					thermos[i].sigmaE()
                )
            );

	Info << "Reading fields Vc\n" << endl;
			VcSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "Vc",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i]
                )
            );

	Info<< "Reading fields Ar\n" << endl;
			ArSolid.set
            (
                i,
                new volVectorField
                (
                    IOobject
                    (
                        "Ar",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i]
                )
            );

	Info<< "Reading fields Ai\n" << endl;
			AiSolid.set
            (
                i,
                new volVectorField
                (
                    IOobject
                    (
                        "Ai",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i]
                )
            );

// 			BSolid.set
// 			(
// 				i,
// 				new volVectorField
// 				(
// 					IOobject
// 					(
// 						"B",
// 						runTime.timeName(),
// 						solidRegions[i],
// 						IOobject::NO_READ,
// 						IOobject::AUTO_WRITE
// 					),
// 					fvc::curl(ArSolid[i] - AiSolid[i])
// 				)
// 			);

// 			EcSolid.set
// 			(
// 				i,
// 				new volVectorField
// 				(
// 					IOobject
// 					(
// 						"Ec",
// 						runTime.timeName(),
// 						solidRegions[i],
// 						IOobject::NO_READ,
// 						IOobject::AUTO_WRITE
// 					),
// 					-fvc::grad(VcSolid[i]) + omegaE*(ArSolid[i] - AiSolid[i])
// 				)
// 			);

// 			JSolid.set
// 			(
// 				i,
// 				new volVectorField
// 				(
// 					IOobject
// 					(
// 						"J",
// 						runTime.timeName(),
// 						solidRegions[i],
// 						IOobject::NO_READ,
// 						IOobject::AUTO_WRITE
// 					),
// 					thermos[i].sigmaE()*EcSolid[i]
// 				)
// 			);

			JcoilSolid.set
			(
				i,
				new volVectorField
				(
					IOobject
					(
						"Jcoil",
						runTime.timeName(),
						solidRegions[i],
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					solidRegions[i]
				)
			);
    }
