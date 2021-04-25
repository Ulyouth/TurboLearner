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
   splineTransportI.H
    splineTransport.C

\*---------------------------------------------------------------------------*/

#ifndef splineTransport_H
#define splineTransport_H

//#include "scalar.H"
#include "FixedList.H"
#include "autoPtr.H"
//#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo>
class splineTransport;

template<class Thermo>
inline splineTransport<Thermo> operator+
(
    const splineTransport<Thermo>&,
    const splineTransport<Thermo>&
);

template<class Thermo>
inline splineTransport<Thermo> operator-
(
    const splineTransport<Thermo>&,
    const splineTransport<Thermo>&
);

template<class Thermo>
inline splineTransport<Thermo> operator*
(
    const scalar,
    const splineTransport<Thermo>&
);

template<class Thermo>
inline splineTransport<Thermo> operator==
(
    const splineTransport<Thermo>&,
    const splineTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const splineTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                        Class icoSpline Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class splineTransport
:
    public Thermo
{
    
        // Public data

        static const int nCoeffs_ = 4;
        typedef FixedList<scalar, nCoeffs_> coeffArray;
        

private:

    // Private data

        // Temperature limits of applicability of functions
        scalar Tlow_, Thigh_, T1_, T2_, T3_, T4_, T5_, T6_, T7_, T8_, T9_, T10_, T11_, T12_, T13_, T14_, T15_;

        coeffArray amuCoeffs_;
        coeffArray bmuCoeffs_;
        coeffArray cmuCoeffs_;
        coeffArray dmuCoeffs_;
        coeffArray emuCoeffs_;
        coeffArray fmuCoeffs_;
        coeffArray gmuCoeffs_;
        coeffArray hmuCoeffs_;
        coeffArray imuCoeffs_;
        coeffArray jmuCoeffs_;
        coeffArray kmuCoeffs_;
        coeffArray lmuCoeffs_;
        coeffArray mmuCoeffs_;
        coeffArray nmuCoeffs_;
        coeffArray omuCoeffs_;
        coeffArray pmuCoeffs_;
        coeffArray akappaCoeffs_;
        coeffArray bkappaCoeffs_;
        coeffArray ckappaCoeffs_;
        coeffArray dkappaCoeffs_;
        coeffArray ekappaCoeffs_;
        coeffArray fkappaCoeffs_;
        coeffArray gkappaCoeffs_;
        coeffArray hkappaCoeffs_;
        coeffArray ikappaCoeffs_;
        coeffArray jkappaCoeffs_;
        coeffArray kkappaCoeffs_;
        coeffArray lkappaCoeffs_;
        coeffArray mkappaCoeffs_;
        coeffArray nkappaCoeffs_;
        coeffArray okappaCoeffs_;
        coeffArray pkappaCoeffs_;        



    // Private Member Functions

        //- Check that input data is valid
        void checkInputData() const;

        //- Return the coefficients corresponding to the given temperature
        inline const coeffArray& mucoeffs(const scalar T) const;

        //- Return the coefficients corresponding to the given temperature
        inline const coeffArray& kappacoeffs(const scalar T) const;

public:

    // Constructors

        //- Construct from components
        inline splineTransport
        (
            const Thermo& t,
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
            const coeffArray& amuCoeffs,
            const coeffArray& bmuCoeffs,
            const coeffArray& cmuCoeffs,
            const coeffArray& dmuCoeffs,
            const coeffArray& emuCoeffs,
            const coeffArray& fmuCoeffs,
            const coeffArray& gmuCoeffs,
            const coeffArray& hmuCoeffs, 
            const coeffArray& imuCoeffs,
            const coeffArray& jmuCoeffs,
            const coeffArray& kmuCoeffs,
            const coeffArray& lmuCoeffs,
            const coeffArray& mmuCoeffs,
            const coeffArray& nmuCoeffs,
            const coeffArray& omuCoeffs,
            const coeffArray& pmuCoeffs, 
            const coeffArray& akappaCoeffs,
            const coeffArray& bkappaCoeffs,
            const coeffArray& ckappaCoeffs,
            const coeffArray& dkappaCoeffs,
            const coeffArray& ekappaCoeffs,
            const coeffArray& fkappaCoeffs,
            const coeffArray& gkappaCoeffs,
            const coeffArray& hkappaCoeffs, 
            const coeffArray& ikappaCoeffs,
            const coeffArray& jkappaCoeffs,
            const coeffArray& kkappaCoeffs,
            const coeffArray& lkappaCoeffs,
            const coeffArray& mkappaCoeffs,
            const coeffArray& nkappaCoeffs,
            const coeffArray& okappaCoeffs,
            const coeffArray& pkappaCoeffs    
           //const Polynomial<PolySize>& rhoPoly
        );

        //- Construct from Istream
        splineTransport(Istream&);

        //- Construct from dictionary
        splineTransport(const dictionary& dict);

        //- Construct as copy
        inline splineTransport(const splineTransport&);

        //- Construct as named copy
        inline splineTransport(const word&, const splineTransport&);

        //- Construct and return a clone
        inline autoPtr<splineTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<splineTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<splineTransport> New(const dictionary& dict);



    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "splineTransport<" + Thermo::typeName() + '>';
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

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& amuCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& bmuCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& cmuCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& dmuCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& emuCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& fmuCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& gmuCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& hmuCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& imuCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& jmuCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& kmuCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& lmuCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& mmuCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& nmuCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& omuCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& pmuCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& akappaCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& bkappaCoeffs() const;
            
            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& ckappaCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& dkappaCoeffs() const;
            
            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& ekappaCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& fkappaCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& gkappaCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& hkappaCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& ikappaCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& jkappaCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& kkappaCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& lkappaCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& mkappaCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& nkappaCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& okappaCoeffs() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& pkappaCoeffs() const;

         //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        //inline splineTransport& operator=(const splineTransport&);
        
        inline void operator+=(const splineTransport&);
        inline void operator-=(const splineTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend splineTransport operator+ <Thermo>
        (
            const splineTransport&,
            const splineTransport&
        );

        friend splineTransport operator- <Thermo>
        (
            const splineTransport&,
            const splineTransport&
        );

        friend splineTransport operator* <Thermo>
        (
            const scalar ,
            const splineTransport&
        );

        friend splineTransport operator== <Thermo>
        (
            const splineTransport&,
            const splineTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const splineTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "splineTransportI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "splineTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
