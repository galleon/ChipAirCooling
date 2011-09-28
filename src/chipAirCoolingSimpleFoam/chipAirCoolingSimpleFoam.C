/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 H. Jasak & H. Rusche
     \\/     M anipulation  | All rights reserved
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    chipAirCoolingFoam

Description
    Steady-State solver for incompressible, turbulent flow of Newtonian fluids
    with conjugate heat transfer, complex heat conduction and radiation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "coupledFvMatrices.H"
#include "regionCouplePolyPatch.H"
#include "radiationModel.H"
#include "thermalModel.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createFluidMesh.H"
#   include "createSolidMesh.H"
#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "createSolidFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // Interpolate to the faces and add thermal resistance
    surfaceScalarField kappaEfff("kappaEfff", fvc::interpolate(kappaEff));
    surfaceScalarField ksolidf("ksolidf", fvc::interpolate(ksolid));

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"
#       include "initConvergenceCheck.H"

        // Detach patches
#       include "detachPatches.H"

        p_rgh.storePrevIter();

#       include "UEqn.H"
#       include "pEqn.H"

        // Update turbulent quantities
        turbulence->correct();

        radiation->correct();

        // Update thermal conductivity in the fluid
        kappaEff = rho*Cp*(turbulence->nu()/Pr + turbulence->nut()/Prt);

        // Update thermal conductivity in the solid
        solidThermo.correct();
        ksolid = solidThermo.k();

        // Coupled patches
#       include "attachPatches.H"

        kappaEff.correctBoundaryConditions();
        ksolid.correctBoundaryConditions();

        //scalarField k1Save = kappaEff.boundaryField()[2];
        //scalarField k2Save = ksolid.boundaryField()[1];

        // Interpolate to the faces and add thermal resistance
        kappaEfff = fvc::interpolate(kappaEff);
        ksolidf = fvc::interpolate(ksolid);
        solidThermo.modifyResistance(ksolidf);


        //kappaEff.boundaryField()[2] = k1Save;
        //ksolid.boundaryField()[1] = k2Save;

        /*
        Info << "k1 = "
            << (scalarField) kappaEfff.boundaryField()[2] << endl;
        Info << "k2 = "
            << (scalarField) ksolidf.boundaryField()[1] << endl;
        Info << "T1 = "
            << (scalarField) T.boundaryField()[2].patchInternalField() << endl;
        Info << "T2 = "
            << (scalarField) Tsolid.boundaryField()[1].patchInternalField() << endl;
        Info << "Tw1 = "
            << (scalarField) T.boundaryField()[2] << endl;
        Info << "Tw2 = "
            << (scalarField) Tsolid.boundaryField()[1] << endl;
            */

#       include "solveEnergy.H"

        /*
        Info << "k1 = "
            << (scalarField) kappaEfff.boundaryField()[2] << endl;
        Info << "k2 = "
            << (scalarField) ksolidf.boundaryField()[1] << endl;
        Info << "T1 = "
            << (scalarField) T.boundaryField()[2].patchInternalField() << endl;
        Info << "T2 = "
            << (scalarField) Tsolid.boundaryField()[1].patchInternalField() << endl;
        Info << "Tw1 = "
            << (scalarField) T.boundaryField()[2] << endl;
        Info << "Tw2 = "
            << (scalarField) Tsolid.boundaryField()[1] << endl;
            */

        volScalarField Tres =
        (
            "Tres",
            rho*Cp*
            (
                fvc::div(phi, T)
            )
            - fvc::laplacian(kappaEfff, T, "laplacian(kappaEff,T)")
        );
        //Tres.writeOpt() = IOobject::AUTO_WRITE;

        volScalarField TsRes =
        (
            "TsRes",
            - fvc::laplacian(ksolidf, Tsolid, "laplacian(k,T)")
            - solidThermo.S()
        );
        //TsRes.writeOpt() = IOobject::AUTO_WRITE;

        Info << "Tres = " << fvc::domainIntegrate(Tres).value() << endl;
        Info << "TsRes = " << fvc::domainIntegrate(TsRes).value() << endl;

        // Update density according to Boussinesq approximation
        rhok = 1.0 - beta*(T - TRef);

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
