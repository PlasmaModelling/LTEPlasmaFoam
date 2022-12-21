/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011, 2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "turbulentPotentialCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentPotentialCoupledMixedFvPatchScalarField::
turbulentPotentialCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    potentialCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    VcnbrName_("undefined-Vcnbr"),
    thicknessLayers_(0),
    sigmaELayers_(0),
    contactRes_(0) //,
//     thermalInertia_(false)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentPotentialCoupledMixedFvPatchScalarField::
turbulentPotentialCoupledMixedFvPatchScalarField
(
    const turbulentPotentialCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    potentialCoupledBase(patch(), psf),
    VcnbrName_(psf.VcnbrName_),
    thicknessLayers_(psf.thicknessLayers_),
    sigmaELayers_(psf.sigmaELayers_),
    contactRes_(psf.contactRes_) //,
//     thermalInertia_(psf.thermalInertia_)
{}


turbulentPotentialCoupledMixedFvPatchScalarField::
turbulentPotentialCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    potentialCoupledBase(patch(), dict),
    VcnbrName_(dict.lookupOrDefault<word>("Vcnbr", "Vc")),
    thicknessLayers_(0),
    sigmaELayers_(0),
    contactRes_(0.0) //,
//     thermalInertia_(dict.lookupOrDefault<Switch>("thermalInertia", false))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
    {
        dict.readEntry("sigmaELayers", sigmaELayers_);

        if (thicknessLayers_.size() > 0)
        {
            // Calculate effective thermal resistance by harmonic averaging
            forAll(thicknessLayers_, iLayer)
            {
                contactRes_ += thicknessLayers_[iLayer]/sigmaELayers_[iLayer];
            }
            contactRes_ = 1.0/contactRes_;
        }
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


turbulentPotentialCoupledMixedFvPatchScalarField::
turbulentPotentialCoupledMixedFvPatchScalarField
(
    const turbulentPotentialCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    potentialCoupledBase(patch(), psf),
    VcnbrName_(psf.VcnbrName_),
    thicknessLayers_(psf.thicknessLayers_),
    sigmaELayers_(psf.sigmaELayers_),
    contactRes_(psf.contactRes_) //,
//     thermalInertia_(psf.thermalInertia_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentPotentialCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const polyMesh& mesh = patch().boundaryMesh().mesh();

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const label patchi = patch().index();
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];


    scalarField Vcc(patchInternalField());
    scalarField& Vcp = *this;

    const turbulentPotentialCoupledMixedFvPatchScalarField&
        nbrField = refCast
            <const turbulentPotentialCoupledMixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>(VcnbrName_)
            );

    // Swap to obtain full local values of neighbour internal field
    scalarField VccNbr(nbrField.patchInternalField());
    mpp.distribute(VccNbr);

    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr;
    if (contactRes_ == 0.0)
    {
        KDeltaNbr = nbrField.sigmaE(nbrField)*nbrPatch.deltaCoeffs();
    }
    else
    {
        KDeltaNbr.setSize(nbrField.size(), contactRes_);
    }
    mpp.distribute(KDeltaNbr);

	scalarField KDelta(sigmaE(Vcp)*patch().deltaCoeffs());

    scalarField qr(Vcp.size(), Zero);

    scalarField qrNbr(Vcp.size(), Zero);

	valueFraction() = KDeltaNbr/(KDeltaNbr + KDelta);
	refValue() = VccNbr;
	refGrad() = (qr + qrNbr)/sigmaE(Vcp);
    mixedFvPatchScalarField::updateCoeffs();

// Info << "Vc values:	" << Vcp << endl;
// Info << "sigmaE values:	" << sigmaE(Vcp) << endl;

    // Restore tag
    UPstream::msgType() = oldTag;
}


void turbulentPotentialCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
	os.writeEntry("Vcnbr", VcnbrName_);
    thicknessLayers_.writeEntry("thicknessLayers", os);
    sigmaELayers_.writeEntry("sigmaELayers", os);

    potentialCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentPotentialCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
