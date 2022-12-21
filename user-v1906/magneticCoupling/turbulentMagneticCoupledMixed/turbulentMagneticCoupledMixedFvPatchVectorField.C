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

#include "turbulentMagneticCoupledMixedFvPatchVectorField.H"
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

turbulentMagneticCoupledMixedFvPatchVectorField::
turbulentMagneticCoupledMixedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    magneticCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    AnbrName_("undefined-Anbr"),
    thicknessLayers_(0),
    sigmaELayers_(0),
    contactRes_(0) //,
//     thermalInertia_(false)
{
    this->refValue() = Zero; //0.0;
    this->refGrad() = Zero; //0.0;
    this->valueFraction() = 1.0;
}


turbulentMagneticCoupledMixedFvPatchVectorField::
turbulentMagneticCoupledMixedFvPatchVectorField
(
    const turbulentMagneticCoupledMixedFvPatchVectorField& psf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(psf, p, iF, mapper),
    magneticCoupledBase(patch(), psf),
    AnbrName_(psf.AnbrName_),
    thicknessLayers_(psf.thicknessLayers_),
    sigmaELayers_(psf.sigmaELayers_),
    contactRes_(psf.contactRes_) //,
//     thermalInertia_(psf.thermalInertia_)
{}


turbulentMagneticCoupledMixedFvPatchVectorField::
turbulentMagneticCoupledMixedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    magneticCoupledBase(patch(), dict),
    AnbrName_(dict.lookupOrDefault<word>("Anbr", "A")),
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

    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = vectorField("refValue", dict, p.size());
        refGrad() = vectorField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = Zero; //0.0;
        valueFraction() = 1.0;
    }
}


turbulentMagneticCoupledMixedFvPatchVectorField::
turbulentMagneticCoupledMixedFvPatchVectorField
(
    const turbulentMagneticCoupledMixedFvPatchVectorField& psf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(psf, iF),
    magneticCoupledBase(patch(), psf),
    AnbrName_(psf.AnbrName_),
    thicknessLayers_(psf.thicknessLayers_),
    sigmaELayers_(psf.sigmaELayers_),
    contactRes_(psf.contactRes_) //,
//     thermalInertia_(psf.thermalInertia_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentMagneticCoupledMixedFvPatchVectorField::updateCoeffs()
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


    vectorField Ac(patchInternalField());
    vectorField& Ap = *this;

    const turbulentMagneticCoupledMixedFvPatchVectorField&
        nbrField = refCast
            <const turbulentMagneticCoupledMixedFvPatchVectorField>
            (
                nbrPatch.lookupPatchField<volVectorField, vector>(AnbrName_)
            );

    // Swap to obtain full local values of neighbour internal field
    vectorField AcNbr(nbrField.patchInternalField());
    mpp.distribute(AcNbr);

    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr;
    if (contactRes_ == 0.0)
    {
        KDeltaNbr = nbrPatch.deltaCoeffs();
    }
    else
    {
        KDeltaNbr.setSize(nbrField.size(), contactRes_);
    }
    mpp.distribute(KDeltaNbr);

	scalarField KDelta(patch().deltaCoeffs());

    vectorField qr(Ap.size(), Zero);

    vectorField qrNbr(Ap.size(), Zero);

	valueFraction() = KDeltaNbr/(KDeltaNbr + KDelta);
	refValue() = AcNbr;
	refGrad() = (qr + qrNbr);
    mixedFvPatchVectorField::updateCoeffs();

// Info << "A values:	" << Ap << endl;
// Info << "sigmaE values:	" << sigmaE(Ap) << endl;

    // Restore tag
    UPstream::msgType() = oldTag;
}


void turbulentMagneticCoupledMixedFvPatchVectorField::write
(
    Ostream& os
) const
{
    mixedFvPatchVectorField::write(os);
	os.writeEntry("Anbr", AnbrName_);
    thicknessLayers_.writeEntry("thicknessLayers", os);
    sigmaELayers_.writeEntry("sigmaELayers", os);

    magneticCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    turbulentMagneticCoupledMixedFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
