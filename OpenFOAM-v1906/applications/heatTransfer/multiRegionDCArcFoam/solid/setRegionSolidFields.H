fvMesh& mesh = solidRegions[i];
solidThermo& thermo = thermos[i];

tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();

tmp<volScalarField> tcp = thermo.Cp();
const volScalarField& cp = tcp();

tmp<volScalarField> tsigmaE = thermo.sigmaE();
const volScalarField& sigmaE = tsigmaE();

const dimensionedScalar& mu0_ = constant::electromagnetic::mu0;

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
volVectorField& A = ASolid[i];
volVectorField& B = BSolid[i];
volVectorField& Ec = EcSolid[i];
volVectorField& J = JSolid[i];

const volScalarField& betav = betavSolid[i];

fv::options& fvOptions = solidHeatSources[i];
