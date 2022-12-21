    const fvMesh& mesh = fluidRegions[i];

    rhoThermo& thermo = thermoFluid[i];
    thermo.validate(args.executable(), "h", "e");

    volScalarField& rho = rhoFluid[i];
	volScalarField& mu = muFluid[i];
	volScalarField& alpha = alphaFluid[i];
	volScalarField& sigmaE = sigmaEFluid[i];
	volScalarField& Cp = CpFluid[i];
	volScalarField& Vc = VcFluid[i];
	volVectorField& A = AFluid[i];
	volScalarField& Jo = JoFluid[i];
	volScalarField& EJ = EJFluid[i];
	volVectorField& Lo = LoFluid[i];
	volVectorField& Ec = EcFluid[i];
	volVectorField& B = BFluid[i];
	volVectorField& J = JFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];

	const dimensionedScalar& mu0_ = constant::electromagnetic::mu0;
	const dimensionedScalar& e_ = constant::electromagnetic::e;
	const dimensionedScalar& kb_ = constant::physicoChemical::k;

    compressible::turbulenceModel& turb = turbulence[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    volScalarField& p_rgh = p_rghFluid[i];

    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    radiation::radiationModel& rad = radiation[i];

    IOMRFZoneList& MRF = MRFfluid[i];
    fv::options& fvOptions = fluidFvOptions[i];

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    bool frozenFlow = frozenFlowFluid[i];

    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];
