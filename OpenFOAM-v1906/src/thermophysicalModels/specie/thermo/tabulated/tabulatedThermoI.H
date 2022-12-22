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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::tabulatedThermo<EquationOfState>::tabulatedThermo
(
    const EquationOfState& st,
	const scalar Hf,
	const scalar Sf,
    const interpolationTable<scalar> cpTable,
    const interpolationTable<scalar> hTable,
	const interpolationTable<scalar> sTable,
	const interpolationTable<scalar> dCpdTTable,
	const interpolationTable<scalar> dHdTTable,
	const interpolationTable<scalar> dSdTTable
// 	const scalar sTab
)
:
    EquationOfState(st),
    Hf_(Hf),
    Sf_(Sf),
    cpTable_(cpTable),
    hTable_(hTable),
    sTable_(sTable),
    dCpdTTable_(dCpdTTable_),
    dHdTTable_(dHdTTable_),
    dSdTTable_(dSdTTable_)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::tabulatedThermo<EquationOfState>::tabulatedThermo
(
    const word& name,
    const tabulatedThermo& ct
)
:
    EquationOfState(name, ct),
    Hf_(ct.Hf_),
    Sf_(ct.Sf_),
    cpTable_(ct.cpTable_),
    hTable_(ct.hTable_),
    sTable_(ct.sTable_),
    dCpdTTable_(ct.dCpdTTable_),
    dHdTTable_(ct.dHdTTable_),
    dSdTTable_(ct.dSdTTable_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::tabulatedThermo<EquationOfState> >
Foam::tabulatedThermo<EquationOfState>::clone() const
{
    return autoPtr<tabulatedThermo<EquationOfState> >
    (
        new tabulatedThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::tabulatedThermo<EquationOfState> >
Foam::tabulatedThermo<EquationOfState>::New(Istream& is)
{
    return autoPtr<tabulatedThermo<EquationOfState> >
    (
        new tabulatedThermo<EquationOfState>(is)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::tabulatedThermo<EquationOfState> >
Foam::tabulatedThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<tabulatedThermo<EquationOfState> >
    (
        new tabulatedThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
//		return cpTable_(p,T, false)*this->W();
		return cpTable_(T) + EquationOfState::Cp(p, T);
}

template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::Ha
(
    const scalar p, const scalar T
) const
{
//	return hTable_(p,T, false)*this->W();
	return hTable_(T) + Hf_ + EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::Hs
(
    const scalar p, const scalar T
) const
{
//	return hTable_(p,T, false)*this->W();
	return Ha(p, T) - Hc();
}

template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::Hc() const
{
	return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::S
(
    const scalar p, const scalar T
) const
{
//     notImplemented
//     (
//         "scalar tabulatedThermo<EquationOfState>::"
//         "S(const scalar p, const scalar T) const"
//     );
//	return T;

	return sTable_(T) + Sf_;
// 	NOTE da rivedere
}

template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::dGdT
(
	const scalar p, const scalar T
) const
{
	return 
	(
		(dHdTTable_(T))
		-T*(dSdTTable_(T))
		-(sTable_(T) + Sf_)
	); //NOTE da rivedere
}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::dCpdT
(
	const scalar p, const scalar T
) const
{
	return dCpdTTable_(T);
}
