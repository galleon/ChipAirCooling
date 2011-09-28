{
    // Solid side
#   include "readSolidControls.H"

    for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        coupledFvScalarMatrix TEqns(2);

        fvScalarMatrix* TFluid = new fvScalarMatrix
        (
            rho*Cp*
            (
                fvm::ddt(T)
              + fvm::div(phi, T)
            )
          - fvm::laplacian(kappaEff, T)
            ==
            radiation->Ru()
          - fvm::Sp(radiation->Rp()*pow3(T), T)
        );

        TFluid->boundaryManipulate(T.boundaryField());
        TFluid->relax(mesh.relaxationFactor("T+T"));

        fvScalarMatrix* TSolid = new fvScalarMatrix
        (
            fvm::ddt(rhoCpsolid, Tsolid)
          - fvm::laplacian(ksolidf, Tsolid, "laplacian(k,T)")
          + fvm::SuSp(-solidThermo.S()/Tsolid, Tsolid)
        );

        TSolid->relax(mesh.relaxationFactor("T+T"));


        // Add fluid equation
        TEqns.set(0, TFluid);

        // Add solid equation
        TEqns.set(1, TSolid);

        TEqns.solve();
    }

    radiation->correct();
}
