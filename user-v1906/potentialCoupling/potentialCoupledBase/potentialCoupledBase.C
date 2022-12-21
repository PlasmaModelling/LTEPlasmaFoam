/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010, 2017-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "potentialCoupledBase.H"
#include "volFields.H"
#include "fluidThermo.H"
#include "solidThermo.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::potentialCoupledBase::SMethodType
>
Foam::potentialCoupledBase::SMethodTypeNames_
{
    { SMethodType::mtFluidThermo, "fluidThermo" },
    { SMethodType::mtSolidThermo, "solidThermo" },
    { SMethodType::mtDirectionalSolidThermo, "directionalSolidThermo" },
    { SMethodType::mtLookup, "lookup" }
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::potentialCoupledBase::potentialCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const word& sigmaEName,
    const word& sigmaEAniName
)
:
    patch_(patch),
    method_(SMethodTypeNames_[calculationType]),
    sigmaEName_(sigmaEName),
    sigmaEAniName_(sigmaEAniName)
{}


Foam::potentialCoupledBase::potentialCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(SMethodTypeNames_.get("sigmaEMethod", dict)),
    sigmaEName_(dict.lookupOrDefault<word>("sigmaE", "none")),
    sigmaEAniName_(dict.lookupOrDefault<word>("sigmaEAni","none"))
{
    switch (method_)
    {
        case mtDirectionalSolidThermo:
        {
            if (!dict.found("sigmaEAni"))
            {
                FatalIOErrorInFunction(dict)
                    << "Did not find entry 'sigmaEAni'"
                       " required for 'sigmaEMethod' "
                    << SMethodTypeNames_[method_]
                    << exit(FatalIOError);
            }

            break;
        }

        case mtLookup:
        {
            if (!dict.found("sigmaE"))
            {
                FatalIOErrorInFunction(dict)
                    << "Did not find entry 'sigmaE'"
                       " required for 'sigmaEMethod' "
                    <<  SMethodTypeNames_[method_] << nl
                    << "    Please set 'sigmaE' to the name of a volScalarField"
                       " or volSymmTensorField"
                    << exit(FatalIOError);
            }

            break;
        }

        default:
        {
            break;
        }
    }
}


Foam::potentialCoupledBase::potentialCoupledBase
(
    const fvPatch& patch,
    const potentialCoupledBase& base
)
:
    patch_(patch),
    method_(base.method_),
    sigmaEName_(base.sigmaEName_),
    sigmaEAniName_(base.sigmaEAniName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::potentialCoupledBase::~potentialCoupledBase()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::potentialCoupledBase::sigmaE
(
    const scalarField& Tp
) const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const label patchi = patch_.index();

    switch (method_)
    {
        case mtFluidThermo:
        {
            typedef compressible::turbulenceModel turbulenceModel;

            const word turbName(turbulenceModel::propertiesName);

            if
            (
                mesh.foundObject<turbulenceModel>(turbName)
            )
            {
                const turbulenceModel& turbModel =
                    mesh.lookupObject<turbulenceModel>(turbName);

//				return turbModel.sigmaE(patchi);
				return turbModel.kappaEff(patchi);
            }
            else if (mesh.foundObject<fluidThermo>(rhoThermo::dictName))
            {
                const fluidThermo& thermo =
                    mesh.lookupObject<fluidThermo>(rhoThermo::dictName);

                return thermo.sigmaE(patchi);
            }
            else if (mesh.foundObject<rhoThermo>(rhoThermo::dictName))
            {
                const rhoThermo& thermo =
                    mesh.lookupObject<rhoThermo>(rhoThermo::dictName);

                return thermo.sigmaE(patchi);
            }
            else if (mesh.foundObject<rhoThermo>("phaseProperties"))
            {
                const rhoThermo& thermo =
                    mesh.lookupObject<rhoThermo>("phaseProperties");

                return thermo.sigmaE(patchi);
            }
            else
            {
                FatalErrorInFunction
                    << "sigmaEMethod defined to employ "
                    << SMethodTypeNames_[method_]
                    << " method, but thermo package not available"
                    << exit(FatalError);
            }

            break;
        }

        case mtSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>(basicThermo::dictName);

            return thermo.sigmaE(patchi);
            break;
        }

        case mtDirectionalSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>(basicThermo::dictName);

            const symmTensorField& sigmaEAni =
                patch_.lookupPatchField<volSymmTensorField, scalar>
                (
                    sigmaEAniName_
                );

            const scalarField& pp = thermo.p().boundaryField()[patchi];

            const symmTensorField sigmaE(sigmaEAni);

            const vectorField n(patch_.nf());

            return n & sigmaE & n;
        }

        case mtLookup:
        {
            if (mesh.foundObject<volScalarField>(sigmaEName_))
            {
                return patch_.lookupPatchField<volScalarField, scalar>
                (
                    sigmaEName_
                );
            }
            else if (mesh.foundObject<volSymmTensorField>(sigmaEName_))
            {
                const symmTensorField& SWall =
                    patch_.lookupPatchField<volSymmTensorField, scalar>
                    (
                        sigmaEName_
                    );

                const vectorField n(patch_.nf());

                return n & SWall & n;
            }
            else
            {
                FatalErrorInFunction
                    << "Did not find field " << sigmaEName_
                    << " on mesh " << mesh.name() << " patch " << patch_.name()
                    << nl
                    << "    Please set 'sigmaE' to the name of a volScalarField"
                    << " or volSymmTensorField."
                    << exit(FatalError);
            }



            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << SMethodTypeNames_[method_] << nl
                << "Please set 'sigmaEMethod' to one of "
                << flatOutput(SMethodTypeNames_.sortedToc()) << nl
                << "and 'sigmaE' to the name of the volScalar"
                << " or volSymmTensor field (if sigmaEMethod=lookup)"
                << exit(FatalError);
        }
    }

    return scalarField(0);
}


void Foam::potentialCoupledBase::write(Ostream& os) const
{
    os.writeEntry("sigmaEMethod", SMethodTypeNames_[method_]);
    os.writeEntry("sigmaE", sigmaEName_);
    os.writeEntry("sigmaEAni", sigmaEAniName_);
}


// ************************************************************************* //
