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
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "chtRcThermalDiffusivityResistanceFvPatchScalarField.H"
#include "chtRcThermalDiffusivitySlaveFvPatchScalarField.H"
#include "chtRcTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

chtRcThermalDiffusivityResistanceFvPatchScalarField::
chtRcThermalDiffusivityResistanceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRegionCoupleBase(p, iF),
    conductivity_(p.size(), 0)
{}


chtRcThermalDiffusivityResistanceFvPatchScalarField::
chtRcThermalDiffusivityResistanceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    chtRegionCoupleBase(p, iF, dict),
    conductivity_("conductivity", dict, p.size())
{}


chtRcThermalDiffusivityResistanceFvPatchScalarField::
chtRcThermalDiffusivityResistanceFvPatchScalarField
(
    const chtRcThermalDiffusivityResistanceFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    chtRegionCoupleBase(ptf, p, iF, mapper),
    conductivity_(ptf.conductivity_)
{}


chtRcThermalDiffusivityResistanceFvPatchScalarField::
chtRcThermalDiffusivityResistanceFvPatchScalarField
(
    const chtRcThermalDiffusivityResistanceFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRegionCoupleBase(ptf, iF),
    conductivity_(ptf.conductivity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void chtRcThermalDiffusivityResistanceFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    calcThermalDiffusivity(*this, shadowPatchField());
}


void
chtRcThermalDiffusivityResistanceFvPatchScalarField::calcThermalDiffusivity
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

    scalarField QrOwn(p.size(), 0.0);
    if (T.radiation())
    {
        QrOwn += p.lookupPatchField<volScalarField, scalar>("Qr");
    }

    scalarField Qr = QrOwn;
    if(T.shadowPatchField().radiation())
    {
        Qr += owner.regionCouplePatch().interpolate
        (
            owner.lookupShadowPatchField<volScalarField, scalar>("Qr")
        );
    }

    scalarField cond(p.size());
    if(isA<chtRcThermalDiffusivitySlaveFvPatchScalarField>(owner))
    {
        cond = owner.regionCouplePatch().interpolate(conductivity_);
    }
    else
    {
        cond = conductivity_;
    }

    //Info << "kOwn = " << kOwn << endl;
    //Info << "kNei = " << kNei << endl;
    //Info << "TcOwn = " << TcOwn << endl;
    //Info << "TcNei = " << TcNei << endl;

    //Info << "QrOwn = " << QrOwn << endl;
    //Info << "Qr = " << Qr << endl;
    //Info << "k1 + k2 = " << (kOwn + kNei) << endl;
 
    //Info << "k = " << k << endl;

    k *= cond/(k*p.deltaCoeffs() + cond);

    //Info << "k = " << k << endl;

    // Correction for radiation and resistance
    scalarField kCorr = kOwn*(Qr + QrOwn*kNei/cond);
    kCorr /= p.deltaCoeffs()*(kOwn + kNei + kOwn*kNei/cond)
        *stabilise(TcNei - TcOwn, SMALL);

    //Info << "kCorr = " << kCorr << endl;

    forAll(k, facei)
    {
        scalar kO = k[facei] + kCorr[facei];

        k[facei] = max(min(kO, 100*k[facei]), 0.01*k[facei]);
    }

    //Info << "k = " << k << endl;
}


void
chtRcThermalDiffusivityResistanceFvPatchScalarField::calcTemperature
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

    scalarField QrOwn(p.size(), 0.0);
    if (owner.radiation())
    {
        QrOwn += p.lookupPatchField<volScalarField, scalar>("Qr");
    }
    
    scalarField Qr = QrOwn;
    if(neighbour.radiation())
    {
        Qr += owner.regionCouplePatch().interpolate
        (
            owner.lookupShadowPatchField<volScalarField, scalar>("Qr")
        );
    }

    scalarField cond(p.size());
    if(isA<chtRcThermalDiffusivitySlaveFvPatchScalarField>(owner))
    {
        cond = owner.regionCouplePatch().interpolate(conductivity_);
    }
    else
    {
        cond = conductivity_;
    }

    scalarField& Tw = owner;

    //Info << "kOwn = " << kOwn << endl;
    //Info << "kNei = " << kNei << endl;
    //Info << "TcOwn = " << TcOwn << endl;
    //Info << "TcNei = " << TcNei << endl;
    //Info << "QrOwn = " << QrOwn << " Sum = " << sum(QrOwn*p.magSf()) << endl;
    //Info << "Qr = " << Qr << " Sum = " << sum(Qr*p.magSf()) << endl;

    Tw = kOwn*TcOwn + kNei*TcNei + Qr + (kOwn*kNei*TcOwn + kNei*QrOwn)/cond;
    Tw /= kOwn + kNei + kOwn*kNei/cond;

    //Info << "Tw = " << Tw << endl;

    //scalarField q1 = (Tw - TcOwn)*kOwn;
    //Info << "q1 = " << q1 << " Sum = " << sum(q1*p.magSf()) << endl;

    //scalarField q2 = (TcNei - TcOwn)*ownerK*p.deltaCoeffs();
    //Info << "q2 = " << q2 << " Sum = " << sum(q2*p.magSf()) << endl;

    owner.fvPatchScalarField::updateCoeffs();
}


void chtRcThermalDiffusivityResistanceFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    conductivity_.writeEntry("conductivity", os);
    os.writeKeyword("remoteField")
        << remoteFieldName() << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    chtRcThermalDiffusivityResistanceFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
