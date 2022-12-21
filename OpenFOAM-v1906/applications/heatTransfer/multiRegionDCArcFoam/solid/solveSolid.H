{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
    {
        fvScalarMatrix hEqn
        (
            (
                thermo.isotropic()
              ? -fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
              : -fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
            )
          ==
            fvOptions(rho, h)
        );

        hEqn.relax();

        fvOptions.constrain(hEqn);

        hEqn.solve();

        fvOptions.correct(h);

		thermo.correct();

		fvScalarMatrix VcEqn
		(
			fvm::laplacian(sigmaE,Vc)
		);

		VcEqn.relax();
		VcEqn.solve();

		fvVectorMatrix AEqn
		(
			fvm::laplacian(A)
			==
			mu0_*sigmaE*fvc::grad(Vc)
		);

		AEqn.relax();
		AEqn.solve();
    }
    
		B = fvc::curl(A);
		Ec = -fvc::grad(Vc);
		J = sigmaE*Ec;
// 		sigmaS = sigmaE;

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
