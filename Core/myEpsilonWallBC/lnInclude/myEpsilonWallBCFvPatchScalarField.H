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

Class
    Foam::myEpsilonWallBCFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition calculated the epsilon wall boundary condition
    according to the Myong Kasagi model.

	 \f[
            epsilon_{wall} = \nu*\frac{\partial^2 k}{\partial y^2}
         \f]

    where

    \vartable
        epsilon | turblence dissipation field
        k       | turblence kinetic energy field
	nu	| kinematic viscosity
	y	| wall normal direction
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            myEpsilonWallBC;
    }
    \endverbatim

See also
    Foam::fixedInternalValueFvPatchField

SourceFiles
    myEpsilonWallBCFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef myEpsilonWallBCFvPatchScalarField_H
#define myEpsilonWallBCFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
           Class myEpsilonWallBCFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class myEpsilonWallBCFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected data

        //- Master patch ID
        label master_;

    // Protected Member Functions

        //- Check the type of the patch
        virtual void checkType();

        //- Set the master patch - master is responsible for updating all
        //  wall function patches
        virtual void setMaster();

        //- Helper function to return non-const access to an epsilon patch
        virtual myEpsilonWallBCFvPatchScalarField& epsilonPatch
        (
            const label patchi
        );

        //- Main driver to calculate epsilon boundary condition
        virtual void calculateEpsilonBoundaryCondition();

        //- Return non-const access to the master patch ID
        virtual label& master()
        {
            return master_;
        }


public:

    //- Runtime type information
    TypeName("myEpsilonWallBC");


    // Constructors

        //- Construct from patch and internal field
        myEpsilonWallBCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myEpsilonWallBCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  myEpsilonWallBCFvPatchScalarField
        //  onto a new patch
        myEpsilonWallBCFvPatchScalarField
        (
            const myEpsilonWallBCFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        myEpsilonWallBCFvPatchScalarField
        (
            const myEpsilonWallBCFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new myEpsilonWallBCFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myEpsilonWallBCFvPatchScalarField
        (
            const myEpsilonWallBCFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new myEpsilonWallBCFvPatchScalarField(*this, iF)
            );
        }

    //- Destructor
    virtual ~myEpsilonWallBCFvPatchScalarField()
    {}


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
