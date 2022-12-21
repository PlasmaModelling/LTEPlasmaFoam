/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

Class
    Foam::magneticCoupledBase

Description
    Common functions used in temperature coupled boundaries.

    The electric conductivity \c sigmaE may be obtained by the following methods:
      - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
        defined by 'sigmaE'
      - 'fluidThermo' : use fluidThermo and default
        compressible::turbulenceModel to calculate sigmaE
      - 'solidThermo' : use solidThermo sigmaE()
      - 'directionalSolidThermo': uses look up for volSymmTensorField for
        transformed sigmaE vector. Field name definable in 'sigmaEAni',
        named 'sigmaEAni' in solid solver by default
      - 'phaseSystem' : used for multiphase thermos

    \par Keywords provided by this class:
      \table
        Property     | Description                | Required    | Default value
        sigmaEMethod  | electric conductivity method        | yes |
        sigmaE        | Name of electric conductivity field | no  | none
        sigmaEAni     | Name of the non-isotropic sigmaE    | no  | sigmaEAni
      \endtable

Usage
    \verbatim
    nonIsotropicWall
    {
        ...
        sigmaEMethod     directionalSolidThermo;
        sigmaE           none;
        sigmaEAni        sigmaEAni;
        ...
    }
    \endverbatim

SourceFiles
    magneticCoupledBase.C

\*---------------------------------------------------------------------------*/

#ifndef magneticCoupledBase_H
#define magneticCoupledBase_H

#include "vectorField.H"
#include "Enum.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class magneticCoupledBase Declaration
\*---------------------------------------------------------------------------*/

class magneticCoupledBase
{
public:

    // Public enumerations

        //- Type of supplied Kappa
        enum SMethodType
        {
            mtFluidThermo,
            mtSolidThermo,
            mtDirectionalSolidThermo,
            mtLookup
        };


protected:

    // Protected data

        static const Enum<SMethodType> SMethodTypeNames_;

        //- Underlying patch
        const fvPatch& patch_;

        //- How to get S
        const SMethodType method_;

        //- Name of thermal conductivity field (if looked up from database)
        const word sigmaEName_;

        //- Name of the non-Isotropic alpha (default: Anialpha)
        const word sigmaEAniName_;


public:

    // Constructors

        //- Construct from patch and K name
        magneticCoupledBase
        (
            const fvPatch& patch,
            const word& calculationMethod,
            const word& sigmaEName,
            const word& sigmaEAniName
        );

        //- Construct from patch and dictionary
        magneticCoupledBase
        (
            const fvPatch& patch,
            const dictionary& dict
        );

         //- Construct from patch and magneticCoupledBase
        magneticCoupledBase
        (
            const fvPatch& patch,
            const magneticCoupledBase& base
        );


    //- Destructor
    virtual ~magneticCoupledBase();

    // Member functions

        //- Method to obtain sigmaE
        word SMethod() const
        {
            return SMethodTypeNames_[method_];
        }

        //- Name of electric conductivity field
        const word& sigmaEName() const
        {
            return sigmaEName_;
        }

        //- Given patch temperature calculate corresponding sigmaE field
        virtual tmp<scalarField> sigmaE(const vectorField& Tp) const;

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
