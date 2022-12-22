{
		volScalarField& he = thermo.he();

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

		B = fvc::curl(A);
		Ec = -fvc::grad(Vc);
		J = sigmaE*Ec;
		Lo = J^B;
		Lo.relax();
		
		Jo = magSqr(J)/sigmaE;
		Jo.relax();
		Info<< "Joule min/max: " << min(Jo).value() << " " << max(Jo).value()<< endl;

		EJ = fvc::div(5*kb_*J*he/(2*e_*Cp));
		EJ.relax();
		Info<< "Electron transport min/max: " << min(mag(EJ)).value() << " " << max(mag(EJ)).value()<< endl;
}
