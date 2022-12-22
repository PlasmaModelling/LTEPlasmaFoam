fvMesh& mesh = solidRegions[i];
solidThermo& thermo = thermos[i];

tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();

tmp<volScalarField> tcp = thermo.Cp();
const volScalarField& cp = tcp();

tmp<volScalarField> tsigmaE = thermo.sigmaE();
const volScalarField& sigmaE = tsigmaE();

const dimensionedScalar& mu0_ = constant::electromagnetic::mu0;

IOdictionary thermophysicalDict
(
    IOobject
    (
        "thermophysicalProperties",    // dictionary name
        runTime.constant(),     // dict is found in "constant"
        mesh,                   // registry for the dict
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

tmp<volSymmTensorField> taniAlpha;
if (!thermo.isotropic())
{
    volSymmTensorField& aniAlpha = aniAlphas[i];
    tmp<volVectorField> tkappaByCp = thermo.Kappa()/cp;
    const coordinateSystem& coodSys = coordinates[i];

    aniAlpha.primitiveFieldRef() =
        coodSys.transformPrincipal
        (
            mesh.cellCentres(),
            tkappaByCp()
        );

    aniAlpha.correctBoundaryConditions();

    taniAlpha = tmp<volSymmTensorField>
    (
        new volSymmTensorField(aniAlpha)
    );
}


volScalarField& h = thermo.he();
volScalarField& Vc = VcSolid[i];
volVectorField& Ar = ArSolid[i];
volVectorField& Ai = AiSolid[i];
// volVectorField& B = BSolid[i];
// volVectorField& Ec = EcSolid[i];
// volVectorField& J = JSolid[i];
volVectorField& Jcoil = JcoilSolid[i];

const volScalarField& betav = betavSolid[i];

fv::options& fvOptions = solidHeatSources[i];
