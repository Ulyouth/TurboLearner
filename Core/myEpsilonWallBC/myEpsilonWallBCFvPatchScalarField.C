/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "myEpsilonWallBCFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H" //for second derivative calculation

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::myEpsilonWallBCFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}

void Foam::myEpsilonWallBCFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<myEpsilonWallBCFvPatchScalarField>(bf[patchi]))
        {
            myEpsilonWallBCFvPatchScalarField& epf = epsilonPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            epf.master() = master;
        }
    }
}

Foam::myEpsilonWallBCFvPatchScalarField&
Foam::myEpsilonWallBCFvPatchScalarField::epsilonPatch(const label patchi)
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const myEpsilonWallBCFvPatchScalarField& epf =
        refCast<const myEpsilonWallBCFvPatchScalarField>(bf[patchi]);

    return const_cast<myEpsilonWallBCFvPatchScalarField&>(epf);
}


void Foam::myEpsilonWallBCFvPatchScalarField::calculateEpsilonBoundaryCondition()
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    DynamicList<label> epsilonPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<myEpsilonWallBCFvPatchScalarField>(bf[patchi]))
        {
            epsilonPatches.append(patchi);
        }
    }

    /* calculate the MK boundary condition for epsilon (epsW = nu*d²k/dy²) */
    const volScalarField& k = db().lookupObject<volScalarField>("k");
    const volScalarField& nu = db().lookupObject<volScalarField>("nu");
    volTensorField secondDrvOfK = (fvc::grad(fvc::grad(k.component(0))));
    volScalarField XXofSecondDrvOfK = secondDrvOfK.component(tensor::XX);
    volScalarField epsilonWall = nu*XXofSecondDrvOfK;

    /* write new boundary condition in boundary field of epsilon */
    forAll(epsilonPatches, i)
    {
	label patchi = epsilonPatches[i];
	myEpsilonWallBCFvPatchScalarField& epf = epsilonPatch(patchi);
	const labelUList& faceCells = bf[patchi].patch().faceCells();
	forAll(faceCells, i)
	{
	    epf[i] = mag(epsilonWall.boundaryField()[patchi][i]); //mag because there exist no negative epsilon values
	}
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myEpsilonWallBCFvPatchScalarField::
myEpsilonWallBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    master_(-1)
{
    checkType();
}


Foam::myEpsilonWallBCFvPatchScalarField::
myEpsilonWallBCFvPatchScalarField
(
    const myEpsilonWallBCFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    master_(-1)
{
    checkType();
}


Foam::myEpsilonWallBCFvPatchScalarField::
myEpsilonWallBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    master_(-1)
{
    checkType();
}


Foam::myEpsilonWallBCFvPatchScalarField::
myEpsilonWallBCFvPatchScalarField
(
    const myEpsilonWallBCFvPatchScalarField& ewfpsf
)
:
    fixedValueFvPatchField<scalar>(ewfpsf),
    master_(-1)
{
    checkType();
}


Foam::myEpsilonWallBCFvPatchScalarField::
myEpsilonWallBCFvPatchScalarField
(
    const myEpsilonWallBCFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ewfpsf, iF),
    master_(-1)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myEpsilonWallBCFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    setMaster();

    if (patch().index() == master_)
    {
	calculateEpsilonBoundaryCondition();
    }

    fvPatchField<scalar>::updateCoeffs();

}

void Foam::myEpsilonWallBCFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        myEpsilonWallBCFvPatchScalarField
    );
}


// ************************************************************************* //
