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
    Foam::icoIAPWS

Description
    Incompressible, polynomial form of equation of state, using a polynomial
    function for density based on IAPWS-IF97.

    WARNING: This is just a dummy class, the functionallity still needs to be
    implemented!

SourceFiles
    icoSplineI.H
    icoSpline.C

\*---------------------------------------------------------------------------*/

#ifndef icoSpline_H
#define icoSpline_H

//#include "scalar.H"
#include "FixedList.H"
#include "autoPtr.H"
//#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie>
class icoSpline;

template<class Specie>
inline icoSpline<Specie> operator+
(
    const icoSpline<Specie>&,
    const icoSpline<Specie>&
);

template<class Specie>
inline icoSpline<Specie> operator-
(
    const icoSpline<Specie>&,
    const icoSpline<Specie>&
);

template<class Specie>
inline icoSpline<Specie> operator*
(
    const scalar,
    const icoSpline<Specie>&
);

template<class Specie>
inline icoSpline<Specie> operator==
(
    const icoSpline<Specie>&,
    const icoSpline<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const icoSpline<Specie>&
);


/*---------------------------------------------------------------------------*\
                        Class icoSpline Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class icoSpline
:
    public Specie
{
    
        // Public data

        static const int nCoeffs_ = 4;
        typedef FixedList<scalar, nCoeffs_> coeffArray;

private:

    // Private data

        // Temperature limits of applicability of functions
        scalar Tlow_, Thigh_, T1_, T2_, T3_, T4_, T5_, T6_, T7_, T8_, T9_, T10_, T11_, T12_, T13_, T14_, T15_, T16_;

        coeffArray arhoCoeffs_;
        coeffArray brhoCoeffs_;
        coeffArray crhoCoeffs_;
        coeffArray drhoCoeffs_;
        coeffArray erhoCoeffs_;
        coeffArray frhoCoeffs_;
        coeffArray grhoCoeffs_;
        coeffArray hrhoCoeffs_;
        coeffArray irhoCoeffs_;
        coeffArray jrhoCoeffs_;
        coeffArray krhoCoeffs_;
        coeffArray lrhoCoeffs_;
        coeffArray mrhoCoeffs_;
        coeffArray nrhoCoeffs_;
        coeffArray orhoCoeffs_;
        coeffArray prhoCoeffs_;
        coeffArray qrhoCoeffs_;



    // Private Member Functions

        //- Check that input data is valid
        void checkInputData() const;

        //- Return the coefficients corresponding to the given temperature
        inline const coeffArray& coeffs(const scalar T) const;

public:

    // Constructors

        //- Construct from components
        inline icoSpline
        (
            const Specie& sp,
            const scalar Tlow,
            const scalar Thigh,
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
            const coeffArray& arhoCoeffs,
            const coeffArray& brhoCoeffs,
            const coeffArray& crhoCoeffs,
            const coeffArray& drhoCoeffs,
            const coeffArray& erhoCoeffs,
            const coeffArray& frhoCoeffs,
            const coeffArray& grhoCoeffs,
            const coeffArray& hrhoCoeffs, 
            const coeffArray& irhoCoeffs,
            const coeffArray& jrhoCoeffs,
            const coeffArray& krhoCoeffs,
            const coeffArray& lrhoCoeffs,
            const coeffArray& mrhoCoeffs,
            const coeffArray& nrhoCoeffs,
            const coeffArray& orhoCoeffs,
            const coeffArray& prhoCoeffs,      
            const coeffArray& qrhoCoeffs 
           //const Polynomial<PolySize>& rhoPoly
        );

        //- Construct from Istream
        icoSpline(Istream&);

        //- Construct from dictionary
        icoSpline(const dictionary& dict);

        //- Construct as copy
        inline icoSpline(const icoSpline&);

        //- Construct as named copy
        inline icoSpline(const word& name, const icoSpline&);

        //- Construct and return a clone
        inline autoPtr<icoSpline> clone() const;

        // Selector from Istream
        inline static autoPtr<icoSpline> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<icoSpline> New(const dictionary& dict);



    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "icoSpline<" + word(Specie::typeName_()) + '>';
        }
        // Access
          
            //- Return const access to the low temperature limit
            inline scalar Tlow() const;

            //- Return const access to the high temperature limit
            inline scalar Thigh() const;

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

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& arhoCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& brhoCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& crhoCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& drhoCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& erhoCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& frhoCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& grhoCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& hrhoCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& irhoCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& jrhoCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& krhoCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& lrhoCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& mrhoCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& nrhoCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& orhoCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& prhoCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& qrhoCoeffs() const;

        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        //inline icoSpline& operator=(const icoSpline&);
        inline void operator+=(const icoSpline&);
        inline void operator-=(const icoSpline&);

        inline void operator*=(const scalar);


    // Friend operators

        friend icoSpline operator+ <Specie>
        (
            const icoSpline&,
            const icoSpline&
        );

        friend icoSpline operator- <Specie>
        (
            const icoSpline&,
            const icoSpline&
        );

        friend icoSpline operator* <Specie>
        (
            const scalar ,
            const icoSpline&
        );

        friend icoSpline operator== <Specie>
        (
            const icoSpline&,
            const icoSpline&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const icoSpline&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "icoSplineI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "icoSpline.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
