/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

Application
    waveFoam

Group
    grpBasicSolvers

Description
	D'Alembert equation solver for the pressure field, under the assumption of uniform speed of sound value.

    \heading Solver details
    The solver is applicable to linear pressure waves diffusion and constant speed of sound coefficient value.  The
    equation is given by:

    \f[
        \left( 1/c^2 ) * \d2d2t{p}  = \div \left( \grad T \right)
    \f]

    Where:
    \vartable
        p     | Pressure field which is solved for
        c   | Speed of sound coefficient
    \endvartable

    \heading Required fields
    \plaintable
    p     | Pressure field which is solved for
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "D'Alembert equation solver for the pressure field, under the assumption of constant speed of sound."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating wave evolution\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
				fvm::d2dt2(p)
			==
				sqr(c)*fvm::laplacian(p)
				
//                 fvm::d2dt2(p) - sqr(c)*fvm::laplacian(p)
//              ==
//                 fvOptions(p)
            );

            fvOptions.constrain(pEqn);
            pEqn.solve();
            fvOptions.correct(p);
        }

        #include "write.H"

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
