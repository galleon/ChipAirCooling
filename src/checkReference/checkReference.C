/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
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

Description
    Sample field data with a choice of interpolation schemes, sampling options
    and write formats.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "fvc.H"

using namespace Foam;

template<class GeoField>
bool checkFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    GeoField* dummy
)
{
    // Search list of objects for volScalarFields
    IOobjectList fieldObjects(objects.lookupClass(GeoField::typeName));

    for
    (
        IOobjectList::iterator iter = fieldObjects.begin();
        iter != fieldObjects.end();
        ++iter
    )
    {
        GeoField current(*iter(), mesh);

        GeoField ref
        (
            IOobject
            (
                "reference"/mesh.time().timeName()/current.name(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        dimensionedScalar small("small", ref.dimensions(), SMALL);

        volScalarField error = mag(current - ref)/stabilise(mag(ref), small);

        dimensionedScalar maxError = max(error);

        Info << ref.name() << tab << maxError.value() << endl;

	if ( maxError.value() > 0.001 )
        {
            return true;
        }
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"

    runTime.setTime(timeDirs[0], 0);

    // Search for list of objects for this time
    IOobjectList objects(mesh, runTime.timeName());

    bool failed =
         checkFields(mesh, objects, reinterpret_cast<volScalarField*>(NULL))
      || checkFields(mesh, objects, reinterpret_cast<volVectorField*>(NULL));

    if ( failed )
    {
        Info << "check failed" << endl;

	return 1;
    }
    else
    {
        Info << "check OK" << endl;

        return 0;
    }
}


// ************************************************************************* //
