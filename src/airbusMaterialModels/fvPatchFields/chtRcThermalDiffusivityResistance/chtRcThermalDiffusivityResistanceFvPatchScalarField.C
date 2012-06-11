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

#include "chtRcThermalDiffusivityResistanceFvPatchScalarField.H"
#include "chtRcThermalDiffusivitySlaveFvPatchScalarField.H"
#include "chtRcTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "harmonic.H"
#include "radiationConstants.H"
#include "VectorN.H"

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

void chtRcThermalDiffusivityResistanceFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchScalarField::evaluate();
}


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
    const fvPatch& p = owner.patch();
    const fvMesh& mesh = p.boundaryMesh().mesh();
    const magLongDelta& mld = magLongDelta::New(mesh);

//    const fvPatch& np = neighbour.patch();
//    const magLongDelta& smld = magLongDelta::New(np.boundaryMesh().mesh());

    const chtRcTemperatureFvPatchScalarField& TwOwn =
        dynamic_cast<const chtRcTemperatureFvPatchScalarField&>
        (
            p.lookupPatchField<volScalarField, scalar>("T")
        );
    scalarField TwNei(p.size(), 0.0);

    scalarField& k = owner;
    const scalarField& fOwn = owner.originalPatchField();
    const scalarField TcOwn = TwOwn.patchInternalField();

    scalarField fNei(p.size());
    scalarField TcNei(p.size());

    scalarField QrOwn(p.size(), 0.0);
    scalarField fourQroOwn(p.size(), 0.0);
    scalarField fourQroNei(p.size());

    if (TwOwn.radiation())
    {
        QrOwn += p.lookupPatchField<volScalarField, scalar>("Qr");
        fourQroOwn += 4.0*radiation::sigmaSB.value()*pow4(TwOwn);
    }

    scalarField Qr = QrOwn;
    scalarField cond(p.size());

    {
        Field<VectorN<scalar, 5> > lData(neighbour.size());

        {
            const scalarField& fNei = neighbour.originalPatchField();
            scalarField TcNei = TwOwn.shadowPatchField().patchInternalField();

            forAll(lData, facei)
            {
                lData[facei][0] = TcNei[facei];
                lData[facei][1] = fNei[facei];
            }
        }

        if(TwOwn.shadowPatchField().radiation())
        {
            const scalarField& TwNei = TwOwn.shadowPatchField();
            const scalarField& QrNei =
                owner.lookupShadowPatchField<volScalarField, scalar>("Qr");

            forAll(lData, facei)
            {
                lData[facei][2] = TwNei[facei];
                lData[facei][3] = QrNei[facei];
            }
        }

        if(isA<chtRcThermalDiffusivitySlaveFvPatchScalarField>(owner))
        {
            forAll(lData, facei)
            {
                lData[facei][4] = conductivity_[facei];
            }
        }

        const Field<VectorN<scalar, 5> > iData =
            owner.regionCouplePatch().interpolate(lData);

        forAll(iData, facei)
        {
            TcNei[facei] = iData[facei][0];
            fNei[facei] = iData[facei][1];
        }

        if(TwOwn.shadowPatchField().radiation())
        {
            forAll(iData, facei)
            {
                TwNei[facei] = iData[facei][2];
                fourQroNei[facei] +=
                    4.0*radiation::sigmaSB.value()*pow4(iData[facei][2]);
                Qr[facei] += iData[facei][3];
            }
        }

        if(isA<chtRcThermalDiffusivitySlaveFvPatchScalarField>(owner))
        {
            forAll(lData, facei)
            {
                cond[facei] = iData[facei][4];
            }
        }
        else
        {
            cond = conductivity_;
        }
    }

    // Do interpolation
    harmonic<scalar> interp(mesh);
    const scalarField weights = interp.weights(fOwn, fNei, p);
    scalarField kHarm = weights*fOwn + (1.0 - weights)*fNei;
    kHarm *= cond/(k*p.deltaCoeffs() + cond);

    const scalarField kOwn = fOwn/(1.0 - p.weights())/mld.magDelta(p.index());
    const scalarField kNei = fNei/p.weights()/mld.magDelta(p.index());

    //Info << "kOwn = " << kOwn << endl;
    //Info << "kNei = " << kNei << endl;
    //Info << "TcOwn = " << TcOwn << endl;
    //Info << "TcNei = " << TcNei << endl;
    //Info << "DeltaT = " << TcNei - TcOwn << endl;

    //Info << "QrOwn = " << QrOwn << endl;
    //Info << "Qr = " << Qr << endl;
    //Info << "kOwn + kNei = " << (kOwn + kNei) << endl;

    scalarField temp = TwNei*(cond + kNei) + fourQroNei;

    k = -TwOwn*
    (
        kOwn*
        (
            TwNei*(cond*(kNei*TcNei + fourQroOwn + fourQroNei + Qr) + kNei*(fourQroOwn + QrOwn))
          + fourQroNei*(fourQroOwn + QrOwn)
        )
        + kOwn*kOwn*TcOwn*temp
    );

    k /= TwOwn*(kOwn*temp + cond*(TwNei*kNei + fourQroNei)) + fourQroOwn*temp;
    k += kOwn*TcOwn;
    k /= stabilise(TcOwn - TcNei, SMALL)*p.deltaCoeffs();

    //Info << "k = " << k << endl;

    forAll(k, facei)
    {
        k[facei] = max(min(k[facei], 100*kHarm[facei]), 0.01*kHarm[facei]);
    }

    //Info << "k = " << k << endl;

    owner.fvPatchScalarField::updateCoeffs();
}


void
chtRcThermalDiffusivityResistanceFvPatchScalarField::calcTemperature
(
    chtRcTemperatureFvPatchScalarField& TwOwn,
    const chtRcTemperatureFvPatchScalarField& TwNei,
    const chtRegionCoupleBase& ownerK
) const
{
    const fvPatch& p = TwOwn.patch();
    const fvMesh& mesh = p.boundaryMesh().mesh();
    const magLongDelta& mld = magLongDelta::New(mesh);

//    const fvPatch& np = neighbour.patch();
//    const magLongDelta& smld = magLongDelta::New(np.boundaryMesh().mesh());

    const scalarField& fOwn = ownerK.originalPatchField();
    const scalarField TcOwn = TwOwn.patchInternalField();

    scalarField fNei(p.size());
    scalarField TcNei(p.size());

    scalarField QrOwn(p.size(), 0.0);
    scalarField fourQroOwn(p.size(), 0.0);
    scalarField fourQroNei(p.size(), 0.0);

    if (TwOwn.radiation())
    {
        QrOwn += p.lookupPatchField<volScalarField, scalar>("Qr");
        fourQroOwn += 4.0*radiation::sigmaSB.value()*pow4(TwOwn);
    }

    scalarField Qr = QrOwn;
    scalarField cond(p.size());

    {
        Field<VectorN<scalar, 5> > lData(TwNei.size());

        {
            const scalarField& fNei =
                ownerK.shadowPatchField().originalPatchField();
            scalarField TcNei = TwOwn.shadowPatchField().patchInternalField();

            forAll(lData, facei)
            {
                lData[facei][0] = TcNei[facei];
                lData[facei][1] = fNei[facei];
            }
        }

        if(TwOwn.shadowPatchField().radiation())
        {
            const scalarField& QrNei =
                TwOwn.lookupShadowPatchField<volScalarField, scalar>("Qr");
            const scalarField& TwNei = TwOwn.shadowPatchField();

            forAll(lData, facei)
            {
                lData[facei][2] = TwNei[facei];
                lData[facei][3] = QrNei[facei];
            }
        }

        if(isA<chtRcThermalDiffusivitySlaveFvPatchScalarField>(ownerK))
        {
            forAll(lData, facei)
            {
                lData[facei][4] = conductivity_[facei];
            }
        }

        const Field<VectorN<scalar, 5> > iData =
            TwOwn.regionCouplePatch().interpolate(lData);

        forAll(iData, facei)
        {
            TcNei[facei] = iData[facei][0];
            fNei[facei] = iData[facei][1];
        }

        if(TwOwn.shadowPatchField().radiation())
        {
            forAll(iData, facei)
            {
                Qr[facei] += iData[facei][3];
                fourQroNei[facei] +=
                    4.0*radiation::sigmaSB.value()*pow4(iData[facei][2]);
            }
        }

        if(isA<chtRcThermalDiffusivitySlaveFvPatchScalarField>(ownerK))
        {
            forAll(lData, facei)
            {
                cond[facei] = iData[facei][4];
            }
        }
        else
        {
            cond = conductivity_;
        }
    }

    // Do interpolation
    harmonic<scalar> interp(mesh);
    scalarField weights = interp.weights(fOwn, fNei, p);
    const scalarField kHarm = weights*fOwn + (1.0 - weights)*fNei;

    const scalarField kOwn = fOwn/(1.0 - p.weights())/mld.magDelta(p.index());
    const scalarField kNei = fNei/p.weights()/mld.magDelta(p.index());

    //Info << "kOwn = " << kOwn << endl;
    //Info << "kNei = " << kNei << endl;
    //Info << "TcOwn = " << TcOwn << endl;
    //Info << "TcNei = " << TcNei << endl;
    //Info << "QrOwn = " << QrOwn << " Sum = " << sum(QrOwn*p.magSf()) << endl;
    //Info << "Qr = " << Qr << " Sum = " << sum(Qr*p.magSf()) << endl;

    scalarField temp = fourQroOwn + QrOwn + kOwn*TcOwn;

    TwOwn *=
    (
        TwNei*(cond*(kNei*TcNei + fourQroOwn + fourQroNei + Qr + kOwn*TcOwn) + kNei*temp)
      + fourQroNei*temp
    )
    /
    (
        TwOwn*((TwNei*(cond*(kOwn + kNei) + kOwn*kNei)) + fourQroNei*(kOwn + cond))
      + fourQroOwn*(TwNei*(cond + kNei) + fourQroNei)
    );

    //Info << "TwOwn = " << TwOwn << endl;

    //scalarField q1 = (TwOwn - TcOwn)*kOwn;
    //Info << "q1 = " << q1 << " Sum = " << sum(q1*p.magSf()) << endl;

    //scalarField q2 = (TcNei - TcOwn)*ownerK*p.deltaCoeffs();
    //Info << "q2 = " << q2 << " Sum = " << sum(q2*p.magSf()) << endl;

    TwOwn.fvPatchScalarField::updateCoeffs();
}


void chtRcThermalDiffusivityResistanceFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    conductivity_.writeEntry("conductivity", os);
    os.writeKeyword("remoteField")
        << remoteFieldName() << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


//- Specify data associated with VectorN<scalar, 5> type is contiguous
template<>
inline bool contiguous<VectorN<scalar, 5> >() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    chtRcThermalDiffusivityResistanceFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
