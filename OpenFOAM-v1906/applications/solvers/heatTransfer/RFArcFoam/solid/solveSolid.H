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

		Jcoil = sigmaE*fvc::grad(Vc);

		fvVectorMatrix ArEqn
		(
			fvm::laplacian(Ar) + mu0_*sigmaE*omegaE*Ai
			==
			- mu0_*Jcoil
		);

		ArEqn.relax();
		ArEqn.solve();

		fvVectorMatrix AiEqn
		(
			fvm::laplacian(Ai) - mu0_*sigmaE*omegaE*Ar
		);

		AiEqn.relax();
		AiEqn.solve();

// 		B = fvc::curl(Ar - Ai);
// 		Ec = -fvc::grad(Vc) + omegaE*(Ar + Ai);
// 		J = sigmaE*Ec;

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
	}
}
