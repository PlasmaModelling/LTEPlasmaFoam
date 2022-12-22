{
// 		volScalarField& he = thermo.he();

// 		fvScalarMatrix VcEqn
// 		(
// 			fvm::laplacian(sigmaE,Vc)
// 		);
// 
// 		VcEqn.relax();
// 		VcEqn.solve();

		fvVectorMatrix ArEqn
		(
			fvm::laplacian(Ar) + mu0_*sigmaE*omegaE*Ai
// 			==
// 			mu0_*sigmaE*fvc::grad(Vc)
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

if (Ecoupling)
{
// 		J = sigmaE*Ec;
// 		Lo = J^B;
		Lo = 0.5*omegaE*sigmaE*(-(Ar^(fvc::curl(Ai))) + (Ai^(fvc::curl(Ar))));
		Lo.relax();

// 		Jo = magSqr(J)/sigmaE;
		Jo = 0.5*sqr(omegaE)*sigmaE*mag(magSqr(Ai) - magSqr(Ar));
		Jo.relax();

		scalar maxx = gMax(Jo);
		if (maxx>=1e+9)
		{
			Jo *= 1e+9/maxx;
		}

		Info<< "Joule min/max: " << min(Jo).value() << " " << max(Jo).value()<< endl;

// 		EJ = fvc::div(5*kb_*J*he/(2*e_*Cp));
// 		EJ.relax();
// 		Info<< "Electron transport min/max: " << min(mag(EJ)).value() << " " << max(mag(EJ)).value()<< endl;
}
else
{
	Info << "Electromagnetic coupling disabled" << endl;
	//NOTE poco elegante ma funzionale. Da cambiare
	Lo *= scalar(0);
	Jo *= scalar(0);
// 	EJ *= scalar(0);
}
}
