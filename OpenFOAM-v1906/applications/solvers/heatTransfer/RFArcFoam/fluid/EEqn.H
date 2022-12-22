{
    volScalarField& he = thermo.he();

	fvScalarMatrix EEqn
	(
		fvm::div(phi, he)
	+ (
			he.name() == "e"
		? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
		: fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
		)
	- fvm::laplacian(turb.alphaEff(), he)
	==
		rho*(U&g)
	+ Jo
// 	+ EJ
	+ rad.Sh(thermo, he)
	+ fvOptions(rho, he)
	);

	EEqn.relax();

	fvOptions.constrain(EEqn);

	EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    rad.correct();

	rho = thermo.rho();
	mu= thermo.mu();
	alpha = thermo.alpha();
	sigmaE= thermo.sigmaE();
	sigmaE.correctBoundaryConditions();
	Cp = thermo.Cp();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
