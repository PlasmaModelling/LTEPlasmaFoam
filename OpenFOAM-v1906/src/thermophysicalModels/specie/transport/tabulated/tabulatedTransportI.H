/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

// #include "interpolation2DTable.H"
#include "interpolationTable.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabulatedTransport<Thermo>::tabulatedTransport
(
    const Thermo& t,
//     const interpolation2DTable<scalar> muTable,
//     const interpolation2DTable<scalar> kappaTable
    const interpolationTable<scalar> muTable,
	const interpolationTable<scalar> kappaTable,
	const interpolationTable<scalar> sigmaETable
)
:
    Thermo(t),
    muTable_(muTable),
    kappaTable_(kappaTable),
    sigmaETable_(sigmaETable)
{}


template<class Thermo>
inline Foam::tabulatedTransport<Thermo>::tabulatedTransport
(
    const word& name,
    const tabulatedTransport& ct
)
:
    Thermo(name, ct),
    muTable_(ct.muTable_),
    kappaTable_(ct.kappaTable_),
    sigmaETable_(ct.sigmaETable_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedTransport<Thermo> >
Foam::tabulatedTransport<Thermo>::clone() const
{
    return autoPtr<tabulatedTransport<Thermo> >
    (
        new tabulatedTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedTransport<Thermo> >
Foam::tabulatedTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<tabulatedTransport<Thermo> >
    (
        new tabulatedTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedTransport<Thermo> >
Foam::tabulatedTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<tabulatedTransport<Thermo> >
    (
        new tabulatedTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::tabulatedTransport<Thermo>::mu
(
    const scalar p, const scalar T
) const
{    
//	return muTable_(p, T, false);
	return muTable_(T);
}


template<class Thermo>
inline Foam::scalar Foam::tabulatedTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
//	return kappaTable_(p, T, false);
	return kappaTable_(T);
}


template<class Thermo>
inline Foam::scalar Foam::tabulatedTransport<Thermo>::alphah
(
    const scalar p, const scalar T
) const
{
//	return kappa(p, T)/this->Cpv(p, T);
// 	return kappa(T)/this->Cpv(T);
	return kappaTable_(T)/this->Cp(p,T);
}

template<class Thermo>
inline Foam::scalar Foam::tabulatedTransport<Thermo>::sigmaE
(
	const scalar p, const scalar T
) const
{
	return sigmaETable_(T);
}
