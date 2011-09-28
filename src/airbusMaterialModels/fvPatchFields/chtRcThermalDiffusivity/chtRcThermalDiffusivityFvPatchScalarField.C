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

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "chtRcThermalDiffusivityFvPatchScalarField.H"
#include "chtRcThermalDiffusivitySlaveFvPatchScalarField.H"
#include "chtRcTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "magLongDelta.H"
#include "radiationConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

chtRcThermalDiffusivityFvPatchScalarField::chtRcThermalDiffusivityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRegionCoupleBase(p, iF)
{}


chtRcThermalDiffusivityFvPatchScalarField::chtRcThermalDiffusivityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    chtRegionCoupleBase(p, iF, dict)
{}


chtRcThermalDiffusivityFvPatchScalarField::chtRcThermalDiffusivityFvPatchScalarField
(
    const chtRcThermalDiffusivityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    chtRegionCoupleBase(ptf, p, iF, mapper)
{}


chtRcThermalDiffusivityFvPatchScalarField::chtRcThermalDiffusivityFvPatchScalarField
(
    const chtRcThermalDiffusivityFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRegionCoupleBase(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
void chtRcThermalDiffusivityFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    Info << "In chtRcThermalDiffusivityFvPatchScalarField::evaluate() on "
        << dimensionedInternalField().name()
        << " in " << patch().boundaryMesh().mesh().name()
        << " " << updated() << endl;

}
    */


void chtRcThermalDiffusivityFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    calcThermalDiffusivity(*this, shadowPatchField());
}


void
chtRcThermalDiffusivityFvPatchScalarField::calcThermalDiffusivity
(
    chtRegionCoupleBase& owner,
    const chtRegionCoupleBase& neighbour
) const
{
    owner.regionCoupleFvPatchScalarField::updateCoeffs();

    const fvPatch& p = owner.patch();
    const magLongDelta& mld = magLongDelta::New(p.boundaryMesh().mesh());

    const fvPatch& np = neighbour.patch();
    const magLongDelta& smld = magLongDelta::New(np.boundaryMesh().mesh());

    const chtRcTemperatureFvPatchScalarField& T =
        dynamic_cast<const chtRcTemperatureFvPatchScalarField&>
        (
            p.lookupPatchField<volScalarField, scalar>("T")
        );

    const scalarField TcOwn = T.patchInternalField();
    const scalarField TcNei = T.patchNeighbourField();

    scalarField& k = owner;
    const scalarField kOwn = owner.originalPatchField()
        /(1 - p.weights())/mld.magDelta(p.index());

    const scalarField kNei = owner.regionCouplePatch().interpolate
    (
        neighbour.originalPatchField()
        /(1 - np.weights())/smld.magDelta(np.index())
    );

    scalarField Qr(p.size(), 0.0);
    if (T.radiation())
    {
        Qr += p.lookupPatchField<volScalarField, scalar>("Qr");
    }

    if(T.shadowPatchField().radiation())
    {
        Qr += owner.regionCouplePatch().interpolate
        (
            owner.lookupShadowPatchField<volScalarField, scalar>("Qr")
        );
    }

    //Info << "kOwn = " << kOwn << endl;
    //Info << "kNei = " << kNei << endl;
    //Info << "TcOwn = " << TcOwn << endl;
    //Info << "TcNei = " << TcNei << endl;
    //Info << "DeltaT = " << TcNei - TcOwn << endl;

    //Info << "QrOwn = " << QrOwn << endl;
    //Info << "QrNei = " << QrNei << endl;
    //Info << "k1 + k2 = " << (kOwn + kNei) << endl;
 
    //Info << "k = " << k << endl;
    
    // Correction for radiation
    //scalarField kCorr = (QrOwn + QrNei)*kOwn;
    //kCorr /= p.deltaCoeffs()*(kOwn + kNei)*stabilise(TcNei - TcOwn, SMALL);

    //scalarField kCorr = kNei*QrOwn - kOwn*QrNei;
    //kCorr /= p.deltaCoeffs()*(kOwn + kNei)*stabilise(TcOwn - TcNei, SMALL);

    scalarField fourQro = 4.0*radiation::sigmaSB.value()*pow4(T);
    fourQro += 4.0*radiation::sigmaSB.value()*pow4
    (
        owner.regionCouplePatch().interpolate(T.shadowPatchField())
    );

    k = T*(fourQro + Qr) - TcOwn*fourQro;
    k /= stabilise(TcNei - TcOwn, SMALL);
    k += T*kNei;
    k *= kOwn;
    k /= p.deltaCoeffs()*(T*(kOwn + kNei) + fourQro);

    //Info << "kCorr = " << kCorr << endl;

    forAll(k, facei)
    {
        const scalarField& kOld = owner.originalPatchField();

        k[facei] = max(min(k[facei], 100*kOld[facei]), 0.01*kOld[facei]);
    }

    //Info << "k = " << k << endl;

    owner.fvPatchScalarField::updateCoeffs();
}


void
chtRcThermalDiffusivityFvPatchScalarField::calcTemperature
(
    chtRcTemperatureFvPatchScalarField& owner,
    const chtRcTemperatureFvPatchScalarField& neighbour,
    const chtRegionCoupleBase& ownerK
) const
{
    const fvPatch& p = owner.patch();
    const magLongDelta& mld = magLongDelta::New(p.boundaryMesh().mesh());

    const fvPatch& np = neighbour.patch();
    const magLongDelta& smld = magLongDelta::New(np.boundaryMesh().mesh());

    const scalarField TcOwn = owner.patchInternalField();
    const scalarField TcNei = owner.patchNeighbourField();

    const scalarField kOwn = ownerK.originalPatchField()
        /(1 - p.weights())/mld.magDelta(p.index());

    const scalarField kNei = owner.regionCouplePatch().interpolate
    (
        ownerK.shadowPatchField().originalPatchField()
        /(1 - np.weights())/smld.magDelta(np.index())
    );

    scalarField Qr(p.size(), 0.0);
    if (owner.radiation())
    {
        Qr += p.lookupPatchField<volScalarField, scalar>("Qr");
    }
    
    if(neighbour.radiation())
    {
        Qr += owner.regionCouplePatch().interpolate
        (
            owner.lookupShadowPatchField<volScalarField, scalar>("Qr")
        );
    }

    scalarField fourQro = 4.0*radiation::sigmaSB.value()*pow4(owner);
    fourQro += 4.0*radiation::sigmaSB.value()*pow4
    (
        owner.regionCouplePatch().interpolate(owner.shadowPatchField())
    );

    scalarField& Tw = owner;

    //Info << "kOwn = " << kOwn << endl;
    //Info << "kNei = " << kNei << endl;
    //Info << "TcOwn = " << TcOwn << endl;
    //Info << "TcNei = " << TcNei << endl;
    //Info << "Qr = " << Qr << " Sum = " << sum(Qr*p.magSf()) << endl;

    //Tw = kOwn*TcOwn + kNei*TcNei + Qr;
    //Tw /= kOwn + kNei;

    // Update is prototional to old Tw! Must divide first
    Tw /= Tw*(kOwn + kNei) + fourQro;
    Tw *= fourQro + Qr + kOwn*TcOwn + kNei*TcNei;

    //Info << "Tw = " << Tw << endl;

    //scalarField q1 = (Tw - TcOwn)*kOwn;
    //Info << "q1 = " << q1 << " Sum = " << sum(q1*p.magSf()) << endl;

    //scalarField q2 = (TcNei - TcOwn)*ownerK*p.deltaCoeffs();
    //Info << "q2 = " << q2 << " Sum = " << sum(q2*p.magSf()) << endl;

    owner.fvPatchScalarField::updateCoeffs();
}


void chtRcThermalDiffusivityFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("remoteField")
        << remoteFieldName() << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    chtRcThermalDiffusivityFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
