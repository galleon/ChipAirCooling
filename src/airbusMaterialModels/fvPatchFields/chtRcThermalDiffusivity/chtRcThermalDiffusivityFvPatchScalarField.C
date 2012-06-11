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
#include "harmonic.H"
#include "radiationConstants.H"
#include "VectorN.H"

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

void chtRcThermalDiffusivityFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchScalarField::evaluate();
}


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

    scalarField& k = owner;
    const scalarField& fOwn = owner.originalPatchField();
    const scalarField TcOwn = TwOwn.patchInternalField();

    scalarField fNei(p.size());
    scalarField TcNei(p.size());

    scalarField Qr(p.size(), 0.0);
    scalarField fourQro(p.size(), 0.0);

    if (TwOwn.radiation())
    {
        Qr += p.lookupPatchField<volScalarField, scalar>("Qr");
        fourQro += 4.0*radiation::sigmaSB.value()*pow4(TwOwn);
    }

    {
        Field<VectorN<scalar, 4> > lData(neighbour.size());

        forAll(lData, facei)
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
            const scalarField& QrNei =
                owner.lookupShadowPatchField<volScalarField, scalar>("Qr");
            const scalarField& TwNei = TwOwn.shadowPatchField();

            forAll(lData, facei)
            {
                lData[facei][2] = TwNei[facei];
                lData[facei][3] = QrNei[facei];
            }
        }

        const Field<VectorN<scalar, 4> > iData =
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
                Qr[facei] += iData[facei][3];
                fourQro[facei] +=
                    4.0*radiation::sigmaSB.value()*pow4(iData[facei][2]);
            }
        }
    }

    // Do interpolation
    harmonic<scalar> interp(mesh);
    const scalarField weights = interp.weights(fOwn, fNei, p);
    const scalarField kHarm = weights*fOwn + (1.0 - weights)*fNei;

    const scalarField kOwn = fOwn/(1.0 - p.weights())/mld.magDelta(p.index());
    const scalarField kNei = fNei/p.weights()/mld.magDelta(p.index());

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

    k = TwOwn*(fourQro + Qr) - TcOwn*fourQro;
    k /= stabilise(TcNei - TcOwn, SMALL);
    k += TwOwn*kNei;
    k *= kOwn;
    k /= p.deltaCoeffs()*(TwOwn*(kOwn + kNei) + fourQro);

    //Info << "k = " << k << endl;

    forAll(k, facei)
    {
        k[facei] = max(min(k[facei], 100*kHarm[facei]), 0.01*kHarm[facei]);
    }

    //Info << "k = " << k << endl;

    owner.fvPatchScalarField::updateCoeffs();
}


void
chtRcThermalDiffusivityFvPatchScalarField::calcTemperature
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

    scalarField Qr(p.size(), 0.0);
    scalarField fourQro(p.size(), 0.0);

    if (TwOwn.radiation())
    {
        Qr += p.lookupPatchField<volScalarField, scalar>("Qr");
        fourQro += 4.0*radiation::sigmaSB.value()*pow4(TwOwn);
    }

    {
        Field<VectorN<scalar, 4> > lData(TwNei.size());

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
            const scalarField& TwNei = TwOwn.shadowPatchField();
            const scalarField& QrNei =
                TwOwn.lookupShadowPatchField<volScalarField, scalar>("Qr");

            forAll(lData, facei)
            {
                lData[facei][2] = TwNei[facei];
                lData[facei][3] = QrNei[facei];
            }
        }

        const Field<VectorN<scalar, 4> > iData =
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
                fourQro[facei] +=
                    4.0*radiation::sigmaSB.value()*pow4(iData[facei][2]);
                Qr[facei] += iData[facei][3];
            }
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
    //Info << "Qr = " << Qr << " Sum = " << sum(Qr*p.magSf()) << endl;

    //TwOwn = kOwn*TcOwn + kNei*TcNei + Qr;
    //TwOwn /= kOwn + kNei;

    // Update is proportional to old Tw! Must divide first
    TwOwn /= TwOwn*(kOwn + kNei) + fourQro;
    TwOwn *= fourQro + Qr + kOwn*TcOwn + kNei*TcNei;

    //Info << "TwOwn = " << TwOwn << endl;

    //scalarField q1 = (TwOwn - TcOwn)*kOwn;
    //Info << "q1 = " << q1 << " Sum = " << sum(q1*p.magSf()) << endl;

    //scalarField q2 = (TcNei - TcOwn)*ownerK*p.deltaCoeffs();
    //Info << "q2 = " << q2 << " Sum = " << sum(q2*p.magSf()) << endl;

    TwOwn.fvPatchScalarField::updateCoeffs();
}


void chtRcThermalDiffusivityFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("remoteField")
        << remoteFieldName() << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


//- Specify data associated with VectorN<scalar, 4> type is contiguous
template<>
inline bool contiguous<VectorN<scalar, 4> >() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    chtRcThermalDiffusivityFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
