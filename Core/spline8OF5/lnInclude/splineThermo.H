/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::janafThermo

Description
    JANAF tables based thermodynamics package templated
    into the equation of state.

SourceFiles
    splineThermoI.H
    splineThermo.C

\*---------------------------------------------------------------------------*/

#ifndef splineThermo_H
#define splineThermo_H

#include "scalar.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class splineThermo;

template<class EquationOfState>
inline splineThermo<EquationOfState> operator+
(
    const splineThermo<EquationOfState>&,
    const splineThermo<EquationOfState>&
);

template<class EquationOfState>
inline splineThermo<EquationOfState> operator-
(
    const splineThermo<EquationOfState>&,
    const splineThermo<EquationOfState>&
);

template<class EquationOfState>
inline splineThermo<EquationOfState> operator*
(
    const scalar,
    const splineThermo<EquationOfState>&
);

template<class EquationOfState>
inline splineThermo<EquationOfState> operator==
(
    const splineThermo<EquationOfState>&,
    const splineThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const splineThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class splineThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class splineThermo
:
    public EquationOfState
{

public:

    // Public data

        static const int nCoeffs_ = 4;
        typedef FixedList<scalar, nCoeffs_> coeffArray;


private:

    // Private data

        // Temperature limits of applicability of functions
        scalar Tlow_, Thigh_,  T0_, T1_, T2_, T3_, T4_, T5_, T6_, T7_, T8_, T9_, T10_, T11_, T12_, T13_, T14_, T15_, T16_, T17_;

        coeffArray a0CpCoeffs_;
        coeffArray aCpCoeffs_;
        coeffArray bCpCoeffs_;
        coeffArray cCpCoeffs_;
        coeffArray dCpCoeffs_;
        coeffArray eCpCoeffs_;
        coeffArray fCpCoeffs_;
        coeffArray gCpCoeffs_;
        coeffArray hCpCoeffs_;
        coeffArray iCpCoeffs_;
        coeffArray jCpCoeffs_;
        coeffArray kCpCoeffs_;
        coeffArray lCpCoeffs_;
        coeffArray mCpCoeffs_;
        coeffArray nCpCoeffs_;
        coeffArray oCpCoeffs_;
        coeffArray pCpCoeffs_;
        coeffArray qCpCoeffs_;
        coeffArray rCpCoeffs_;
        coeffArray a0HsCoeffs_;
        coeffArray aHsCoeffs_;
        coeffArray bHsCoeffs_;
        coeffArray cHsCoeffs_;
        coeffArray dHsCoeffs_;
        coeffArray eHsCoeffs_;
        coeffArray fHsCoeffs_;
        coeffArray gHsCoeffs_;
        coeffArray hHsCoeffs_;
        coeffArray iHsCoeffs_;
        coeffArray jHsCoeffs_;
        coeffArray kHsCoeffs_;
        coeffArray lHsCoeffs_;
        coeffArray mHsCoeffs_;
        coeffArray nHsCoeffs_;
        coeffArray oHsCoeffs_;
        coeffArray pHsCoeffs_;
        coeffArray qHsCoeffs_;
        coeffArray rHsCoeffs_;


    // Private Member Functions

        //- Check that input data is valid
        void checkInputData() const;

        //- Return the coefficients corresponding to the given temperature
        inline const coeffArray& Cpcoeffs(const scalar T) const;

        //- Return the coefficients corresponding to the given temperature
        inline const coeffArray& Hscoeffs(const scalar T) const;


public:

    // Constructors

        //- Construct from components
        inline splineThermo
        (
            const EquationOfState& st,
            const scalar Tlow,
            const scalar Thigh,
            const scalar T0,
            const scalar T1,
            const scalar T2,
            const scalar T3,
            const scalar T4,
            const scalar T5,
            const scalar T6,
            const scalar T7,
            const scalar T8,
            const scalar T9,
            const scalar T10,
            const scalar T11,
            const scalar T12,
            const scalar T13,
            const scalar T14,
            const scalar T15,
            const scalar T16,
            const scalar T17,
            const coeffArray& a0CpCoeffs,
            const coeffArray& aCpCoeffs,
            const coeffArray& bCpCoeffs,
            const coeffArray& cCpCoeffs,
            const coeffArray& dCpCoeffs,
            const coeffArray& eCpCoeffs,
            const coeffArray& fCpCoeffs,
            const coeffArray& gCpCoeffs,
            const coeffArray& hCpCoeffs, 
            const coeffArray& iCpCoeffs,
            const coeffArray& jCpCoeffs,
            const coeffArray& kCpCoeffs,
            const coeffArray& lCpCoeffs, 
            const coeffArray& mCpCoeffs,
            const coeffArray& nCpCoeffs,
            const coeffArray& oCpCoeffs,
            const coeffArray& pCpCoeffs,     
            const coeffArray& qCpCoeffs,
            const coeffArray& rCpCoeffs,
            const coeffArray& a0HsCoeffs,
            const coeffArray& aHsCoeffs,
            const coeffArray& bHsCoeffs,
            const coeffArray& cHsCoeffs,
            const coeffArray& dHsCoeffs,
            const coeffArray& eHsCoeffs,
            const coeffArray& fHsCoeffs,
            const coeffArray& gHsCoeffs,
            const coeffArray& hHsCoeffs, 
            const coeffArray& iHsCoeffs,
            const coeffArray& jHsCoeffs,
            const coeffArray& kHsCoeffs,
            const coeffArray& lHsCoeffs, 
            const coeffArray& mHsCoeffs,
            const coeffArray& nHsCoeffs,
            const coeffArray& oHsCoeffs,
            const coeffArray& pHsCoeffs,     
            const coeffArray& qHsCoeffs,
            const coeffArray& rHsCoeffs
        );

        //- Construct from Istream
        splineThermo(Istream&);

        //- Construct from dictionary
        splineThermo(const dictionary& dict);

        //- Construct as a named copy
        inline splineThermo(const word&, const splineThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "spline<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Access

            //- Return const access to the low temperature limit
            inline scalar Tlow() const;

            //- Return const access to the high temperature limit
            inline scalar Thigh() const;


            //- Return const access to the common temperature
            inline scalar T0() const;


            //- Return const access to the common temperature
            inline scalar T1() const;

            //- Return const access to the common temperature
            inline scalar T2() const;

            //- Return const access to the common temperature
            inline scalar T3() const;

            //- Return const access to the common temperature
            inline scalar T4() const;

            //- Return const access to the common temperature
            inline scalar T5() const;

            //- Return const access to the common temperature
            inline scalar T6() const;

            //- Return const access to the common temperature
            inline scalar T7() const;

            //- Return const access to the common temperature
            inline scalar T8() const;

            //- Return const access to the common temperature
            inline scalar T9() const;

            //- Return const access to the common temperature
            inline scalar T10() const;

            //- Return const access to the common temperature
            inline scalar T11() const;

            //- Return const access to the common temperature
            inline scalar T12() const;

            //- Return const access to the common temperature
            inline scalar T13() const;

            //- Return const access to the common temperature
            inline scalar T14() const;

            //- Return const access to the common temperature
            inline scalar T15() const;

            //- Return const access to the common temperature
            inline scalar T16() const;

            //- Return const access to the common temperature
            inline scalar T17() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& a0CpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& aCpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& bCpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& cCpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& dCpCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& eCpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& fCpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& gCpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& hCpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& iCpCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& jCpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& kCpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& lCpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& mCpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& nCpCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& oCpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& pCpCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& qCpCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& rCpCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& a0HsCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& aHsCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& bHsCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& cHsCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& dHsCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& eHsCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& fHsCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& gHsCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& hHsCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& iHsCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& jHsCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& kHsCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& lHsCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& mHsCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& nHsCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& oHsCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& pHsCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& qHsCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& rHsCoeffs() const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kg K)]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const splineThermo&);
        inline void operator-=(const splineThermo&);


    // Friend operators

        friend splineThermo operator+ <EquationOfState>
        (
            const splineThermo&,
            const splineThermo&
        );

        friend splineThermo operator- <EquationOfState>
        (
            const splineThermo&,
            const splineThermo&
        );

        friend splineThermo operator* <EquationOfState>
        (
            const scalar,
            const splineThermo&
        );

        friend splineThermo operator== <EquationOfState>
        (
            const splineThermo&,
            const splineThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const splineThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "splineThermoI.H"

#ifdef NoRepository
#   include "splineThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
