    if (runTime.writeTime())
    {
        volVectorField gradp(fvc::grad(p));

        volScalarField gradpx
        (
            IOobject
            (
                "gradpx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradp.component(vector::X)
        );

        volScalarField gradpy
        (
            IOobject
            (
                "gradpy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradp.component(vector::Y)
        );

        volScalarField gradpz
        (
            IOobject
            (
                "gradpz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradp.component(vector::Z)
        );


        runTime.write();
    }
