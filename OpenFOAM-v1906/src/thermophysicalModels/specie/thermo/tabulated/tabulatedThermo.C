/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "tabulatedThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::tabulatedThermo<EquationOfState>::tabulatedThermo(Istream& is)
:
    EquationOfState(is),
    cpTable_(readScalar(is)),
    hTable_(readScalar(is)),
    sTable_(readScalar(is)),
    dCpdTTable_(readScalar(is)),
    dHdTTable_(readScalar(is)),
    dSdTTable_(readScalar(is))
{
    is.check("tabulatedThermo::tabulatedThermo(Istream& is)");
}


template<class EquationOfState>
Foam::tabulatedThermo<EquationOfState>::tabulatedThermo(const dictionary& dict)
:
    EquationOfState(dict),
    Hf_(dict.subDict("thermodynamics").get<scalar>("Hf")),
    Sf_(dict.subDict("thermodynamics").get<scalar>("Sf")),
    cpTable_ ("constant/Cp"),
    hTable_("constant/H"),
    sTable_("constant/S"),
    dCpdTTable_("constant/dCpdT"),
    dHdTTable_("constant/dHdT"),
    dSdTTable_("constant/dSdT")
{ //NOTE da sistemare

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::tabulatedThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

//     dictionary dict("thermodynamics");
//     os  << indent << dict.dictName() << dict;
	
	os.beginBlock("thermodynamics");
	os.writeEntry("Hf", Hf_);
	os.writeEntry("Sf", Sf_);
	os.endBlock();
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tabulatedThermo<EquationOfState>& ct
)
{
//     os  << static_cast<const EquationOfState&>(ct) << tab
//         << ct.cpTable << tab << ct.hTable;
	
		os  << static_cast<const EquationOfState&>(ct) << tab
		<< ct.cpTable << tab << ct.hTable << tab << ct.sTable;

    os.check("Ostream& operator<<(Ostream& os, const tabulatedThermo& ct)");
    return os;
}


// ************************************************************************* //
