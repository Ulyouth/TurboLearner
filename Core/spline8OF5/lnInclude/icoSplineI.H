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

\*---------------------------------------------------------------------------*/

#include "icoSpline.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::icoSpline<Specie>::icoSpline
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
    const typename icoSpline<Specie>::coeffArray& arhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& brhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& crhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& drhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& erhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& frhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& grhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& hrhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& irhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& jrhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& krhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& lrhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& mrhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& nrhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& orhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& prhoCoeffs,
    const typename icoSpline<Specie>::coeffArray& qrhoCoeffs
)
:
    Specie(sp),
    Tlow_(Tlow),
    Thigh_(Thigh),
    T1_(T1),
    T2_(T2),
    T3_(T3),
    T4_(T4),
    T5_(T5),
    T6_(T6),
    T7_(T7),
    T8_(T8),
    T9_(T9),
    T10_(T10),
    T11_(T11),
    T12_(T12),
    T13_(T13),
    T14_(T14),
    T15_(T15),
    T16_(T16)

{
    for (register label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        arhoCoeffs_[coefLabel] = arhoCoeffs[coefLabel];
        brhoCoeffs_[coefLabel] = brhoCoeffs[coefLabel];
        crhoCoeffs_[coefLabel] = crhoCoeffs[coefLabel];
        drhoCoeffs_[coefLabel] = drhoCoeffs[coefLabel];
        erhoCoeffs_[coefLabel] = erhoCoeffs[coefLabel];
        frhoCoeffs_[coefLabel] = frhoCoeffs[coefLabel];
        grhoCoeffs_[coefLabel] = grhoCoeffs[coefLabel];
        hrhoCoeffs_[coefLabel] = hrhoCoeffs[coefLabel];
        irhoCoeffs_[coefLabel] = irhoCoeffs[coefLabel];
        jrhoCoeffs_[coefLabel] = jrhoCoeffs[coefLabel];
        krhoCoeffs_[coefLabel] = krhoCoeffs[coefLabel];
        lrhoCoeffs_[coefLabel] = lrhoCoeffs[coefLabel];
        mrhoCoeffs_[coefLabel] = mrhoCoeffs[coefLabel];
        nrhoCoeffs_[coefLabel] = nrhoCoeffs[coefLabel];
        orhoCoeffs_[coefLabel] = orhoCoeffs[coefLabel];
        prhoCoeffs_[coefLabel] = prhoCoeffs[coefLabel];
        qrhoCoeffs_[coefLabel] = qrhoCoeffs[coefLabel];

    }
}


template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::coeffs
(
    const scalar T
) const
{
    if (T < T1_)
    {
        return arhoCoeffs_;
    }
    else if (T < T2_)
    {
        return brhoCoeffs_;
    }
    else if (T < T3_)
    {
        return crhoCoeffs_;
    }    
    else if ( T < T4_)
    {
        return drhoCoeffs_;
    }
    else if (T < T5_)
    {
        return erhoCoeffs_;
    }
    else if (T < T6_)
    {
        return frhoCoeffs_;
    }
    else if (T < T7_)
    {
        return grhoCoeffs_;
    }    
    else if ( T < T8_)
    {
        return hrhoCoeffs_;
    }
    else if (T < T9_)
    {
        return irhoCoeffs_;
    }
    else if (T < T10_)
    {
        return jrhoCoeffs_;
    }
    else if (T < T11_)
    {
        return krhoCoeffs_;
    }    
    else if ( T < T12_)
    {
        return lrhoCoeffs_;
    }
    else if (T < T13_)
    {
        return mrhoCoeffs_;
    }
    else if (T < T14_)
    {
        return nrhoCoeffs_;
    }
    else if (T < T15_)
    {
        return orhoCoeffs_;
    }    
    else if ( T < T16_)
    {
        return prhoCoeffs_;
    }
    else
    {
        return qrhoCoeffs_;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::icoSpline<Specie>::icoSpline
(
    const word& name,
    const icoSpline<Specie>& iS
)
:
    Specie(name, iS),
    Tlow_(iS.Tlow_),
    Thigh_(iS.Thigh_),
    T1_(iS.T1_),
    T2_(iS.T2_),
    T3_(iS.T3_),
    T4_(iS.T4_),
    T5_(iS.T5_),
    T6_(iS.T6_),
    T7_(iS.T7_),
    T8_(iS.T8_),
    T9_(iS.T9_),
    T10_(iS.T10_),
    T11_(iS.T11_),
    T12_(iS.T12_),
    T13_(iS.T13_),
    T14_(iS.T14_),
    T15_(iS.T15_),
    T16_(iS.T16_)
{
    for (register label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        arhoCoeffs_[coefLabel] = iS.arhoCoeffs_[coefLabel];
        brhoCoeffs_[coefLabel] = iS.brhoCoeffs_[coefLabel];
        crhoCoeffs_[coefLabel] = iS.crhoCoeffs_[coefLabel];
        drhoCoeffs_[coefLabel] = iS.drhoCoeffs_[coefLabel];
        erhoCoeffs_[coefLabel] = iS.erhoCoeffs_[coefLabel];
        frhoCoeffs_[coefLabel] = iS.frhoCoeffs_[coefLabel];
        grhoCoeffs_[coefLabel] = iS.grhoCoeffs_[coefLabel];
        hrhoCoeffs_[coefLabel] = iS.hrhoCoeffs_[coefLabel];
        irhoCoeffs_[coefLabel] = iS.irhoCoeffs_[coefLabel];
        jrhoCoeffs_[coefLabel] = iS.jrhoCoeffs_[coefLabel];
        krhoCoeffs_[coefLabel] = iS.krhoCoeffs_[coefLabel];
        lrhoCoeffs_[coefLabel] = iS.lrhoCoeffs_[coefLabel];
        mrhoCoeffs_[coefLabel] = iS.mrhoCoeffs_[coefLabel];
        nrhoCoeffs_[coefLabel] = iS.nrhoCoeffs_[coefLabel];
        orhoCoeffs_[coefLabel] = iS.orhoCoeffs_[coefLabel];
        prhoCoeffs_[coefLabel] = iS.prhoCoeffs_[coefLabel];
        qrhoCoeffs_[coefLabel] = iS.qrhoCoeffs_[coefLabel];
    }
}

//template<class Specie>
//inline Foam::autoPtr<Foam::icoSpline<Specie> >
//Foam::icoSpline<Specie>::clone() const
//{
//    return autoPtr<icoSpline<Specie> >
//    (
//        new icoSpline<Specie>(*this)
//    );//
//}

//template<class Specie>
//inline Foam::autoPtr<Foam::icoSpline<Specie> >
//Foam::icoSpline<Specie>::New(Istream& is)
//{
//    return autoPtr<icoSpline<Specie> >
//    (
//        new icoSpline<Specie>(is)
//    );
//}


//template<class Specie>
//inline Foam::autoPtr<Foam::icoSpline<Specie> >
//Foam::icoSpline<Specie>::New(const dictionary& dict)
//{
//    return autoPtr<icoSpline<Specie> >
//    (
//        new icoSpline<Specie>(dict)
//    );
//}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::Thigh() const
{
    return Thigh_;
}


template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T1() const
{
    return T1_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T2() const
{
    return T2_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T3() const
{
    return T3_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T4() const
{
    return T4_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T5() const
{
    return T5_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T6() const
{
    return T6_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T7() const
{
    return T7_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T8() const
{
    return T8_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T9() const
{
    return T9_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T10() const
{
    return T10_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T11() const
{
    return T11_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T12() const
{
    return T12_;
}
template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T13() const
{
    return T13_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T14() const
{
    return T14_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T15() const
{
    return T15_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::T16() const
{
    return T16_;
}
template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::arhoCoeffs() const
{
    return arhoCoeffs_;
}

template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::brhoCoeffs() const
{
    return brhoCoeffs_;
}

template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::crhoCoeffs() const
{
    return crhoCoeffs_;
}

template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::drhoCoeffs() const
{
    return drhoCoeffs_;
}
template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::erhoCoeffs() const
{
    return erhoCoeffs_;
}
template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::frhoCoeffs() const
{
    return frhoCoeffs_;
}

template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::grhoCoeffs() const
{
    return grhoCoeffs_;
}

template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::hrhoCoeffs() const
{
    return hrhoCoeffs_;
}

template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::irhoCoeffs() const
{
    return irhoCoeffs_;
}
template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::jrhoCoeffs() const
{
    return jrhoCoeffs_;
}
template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::krhoCoeffs() const
{
    return krhoCoeffs_;
}

template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::lrhoCoeffs() const
{
    return lrhoCoeffs_;
}

template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::mrhoCoeffs() const
{
    return mrhoCoeffs_;
}

template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::nrhoCoeffs() const
{
    return nrhoCoeffs_;
}
template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::orhoCoeffs() const
{
    return orhoCoeffs_;
}
template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::prhoCoeffs() const
{
    return prhoCoeffs_;
}

template<class Specie>
inline const typename Foam::icoSpline<Specie>::coeffArray&
Foam::icoSpline<Specie>::qrhoCoeffs() const
{
    return qrhoCoeffs_;
}

template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::rho
(
    scalar p,
    scalar T
) const
{
     const coeffArray& a = coeffs(T);
    return  a[0]*pow((T-273.15),3)+a[1]*pow((T-273.15),2)+a[2]*pow((T-273.15),1)+a[3];
    //RR*
    //(
    //    ((((a[4]/5.0*T + a[3]/4.0)*T + a[2]/3.0)*T + a[1]/2.0)*T + a[0])*T
    //  + a[5]
    //);//
}


template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    return 0.0; 
}


template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    return 0.0; // TODO: needs changing
}


template<class Specie>
inline Foam::scalar Foam::icoSpline<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return 0.0; // TODO: needs changing
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


template<class Specie>
inline void Foam::icoSpline<Specie>::operator+=
(
    const icoSpline<Specie>& iS
)
{
    scalar molr1 = this->nMoles();

    Specie::operator+=(iS);

    molr1 /= this->nMoles();
    scalar molr2 = iS.nMoles()/this->nMoles();

    Tlow_ = max(Tlow_, iS.Tlow_);
    Thigh_ = min(Thigh_, iS.Thigh_);

    if (icoSpline<Specie>::debug && notEqual(T1_, iS.T1_))
    {
        FatalErrorIn
        (
            "icoSpline<Specie>::operator+="
            "(const icoSpline<Specie>& iS) const"
        )   << "T1 " << T1_ << " for "
            << (this->name().size() ? this->name() : "others")
            << " != " << iS.T1_ << " for "
            << (iS.name().size() ? iS.name() : "others")
            << exit(FatalError);
    }

    for
    (
        register label coefLabel=0;
        coefLabel<icoSpline<Specie>::nCoeffs_;
        coefLabel++
    )
    {
        arhoCoeffs_[coefLabel] =
            molr1*arhoCoeffs_[coefLabel]
          + molr2*iS.arhoCoeffs_[coefLabel];

        brhoCoeffs_[coefLabel] =
            molr1*brhoCoeffs_[coefLabel]
          + molr2*iS.brhoCoeffs_[coefLabel];
   
        crhoCoeffs_[coefLabel] =
            molr1*crhoCoeffs_[coefLabel]
          + molr2*iS.crhoCoeffs_[coefLabel];

        drhoCoeffs_[coefLabel] =
            molr1*drhoCoeffs_[coefLabel]
          + molr2*iS.drhoCoeffs_[coefLabel];
     
        erhoCoeffs_[coefLabel] =
            molr1*erhoCoeffs_[coefLabel]
          + molr2*iS.erhoCoeffs_[coefLabel];

        frhoCoeffs_[coefLabel] =
            molr1*frhoCoeffs_[coefLabel]
          + molr2*iS.frhoCoeffs_[coefLabel];

        grhoCoeffs_[coefLabel] =
            molr1*grhoCoeffs_[coefLabel]
          + molr2*iS.grhoCoeffs_[coefLabel];
   
        hrhoCoeffs_[coefLabel] =
            molr1*hrhoCoeffs_[coefLabel]
          + molr2*iS.hrhoCoeffs_[coefLabel];

        irhoCoeffs_[coefLabel] =
            molr1*irhoCoeffs_[coefLabel]
          + molr2*iS.irhoCoeffs_[coefLabel];
     
        jrhoCoeffs_[coefLabel] =
            molr1*jrhoCoeffs_[coefLabel]
          + molr2*iS.jrhoCoeffs_[coefLabel];

        krhoCoeffs_[coefLabel] =
            molr1*krhoCoeffs_[coefLabel]
          + molr2*iS.krhoCoeffs_[coefLabel];

        lrhoCoeffs_[coefLabel] =
            molr1*lrhoCoeffs_[coefLabel]
          + molr2*iS.lrhoCoeffs_[coefLabel];
   
        mrhoCoeffs_[coefLabel] =
            molr1*mrhoCoeffs_[coefLabel]
          + molr2*iS.mrhoCoeffs_[coefLabel];

        nrhoCoeffs_[coefLabel] =
            molr1*nrhoCoeffs_[coefLabel]
          + molr2*iS.nrhoCoeffs_[coefLabel];
     
        orhoCoeffs_[coefLabel] =
            molr1*orhoCoeffs_[coefLabel]
          + molr2*iS.orhoCoeffs_[coefLabel];

        prhoCoeffs_[coefLabel] =
            molr1*prhoCoeffs_[coefLabel]
          + molr2*iS.prhoCoeffs_[coefLabel];

        qrhoCoeffs_[coefLabel] =
            molr1*qrhoCoeffs_[coefLabel]
          + molr2*iS.qrhoCoeffs_[coefLabel];
   
    }
}



template<class Specie>
inline void Foam::icoSpline<Specie>::operator-=
(
    const icoSpline<Specie>& iS
)
{
    scalar molr1 = this->nMoles();

    Specie::operator-=(iS);

    molr1 /= this->nMoles();
    scalar molr2 = iS.nMoles()/this->nMoles();

    Tlow_ = max(Tlow_, iS.Tlow_);
    Thigh_ = min(Thigh_, iS.Thigh_);

    if (icoSpline<Specie>::debug && notEqual(T1_, iS.T1_))
    {
        FatalErrorIn
        (
            "icoSpline<Specie>::operator+="
            "(const icoSpline<Specie>& is) const"
        )   << "T1 " << T1_ << " for "
            << (this->name().size() ? this->name() : "others")
            << " != " << iS.T1_ << " for "
            << (iS.name().size() ? iS.name() : "others")
            << exit(FatalError);
    }

    for
    (
        register label coefLabel=0;
        coefLabel<icoSpline<Specie>::nCoeffs_;
        coefLabel++
    )
    {
        arhoCoeffs_[coefLabel] =
            molr1*arhoCoeffs_[coefLabel]
          - molr2*iS.arhoCoeffs_[coefLabel];

        brhoCoeffs_[coefLabel] =
            molr1*brhoCoeffs_[coefLabel]
          - molr2*iS.brhoCoeffs_[coefLabel];
   
        crhoCoeffs_[coefLabel] =
            molr1*crhoCoeffs_[coefLabel]
          - molr2*iS.crhoCoeffs_[coefLabel];

        drhoCoeffs_[coefLabel] =
            molr1*drhoCoeffs_[coefLabel]
          - molr2*iS.drhoCoeffs_[coefLabel];
     
        erhoCoeffs_[coefLabel] =
            molr1*erhoCoeffs_[coefLabel]
          - molr2*iS.erhoCoeffs_[coefLabel];

        frhoCoeffs_[coefLabel] =
            molr1*frhoCoeffs_[coefLabel]
          - molr2*iS.frhoCoeffs_[coefLabel];

        grhoCoeffs_[coefLabel] =
            molr1*grhoCoeffs_[coefLabel]
          - molr2*iS.grhoCoeffs_[coefLabel];
   
        hrhoCoeffs_[coefLabel] =
            molr1*hrhoCoeffs_[coefLabel]
          - molr2*iS.hrhoCoeffs_[coefLabel];

        irhoCoeffs_[coefLabel] =
            molr1*irhoCoeffs_[coefLabel]
          - molr2*iS.irhoCoeffs_[coefLabel];
     
        jrhoCoeffs_[coefLabel] =
            molr1*jrhoCoeffs_[coefLabel]
          - molr2*iS.jrhoCoeffs_[coefLabel];

        krhoCoeffs_[coefLabel] =
            molr1*krhoCoeffs_[coefLabel]
          - molr2*iS.krhoCoeffs_[coefLabel];

        lrhoCoeffs_[coefLabel] =
            molr1*lrhoCoeffs_[coefLabel]
          - molr2*iS.lrhoCoeffs_[coefLabel];
   
        mrhoCoeffs_[coefLabel] =
            molr1*mrhoCoeffs_[coefLabel]
          - molr2*iS.mrhoCoeffs_[coefLabel];

        nrhoCoeffs_[coefLabel] =
            molr1*nrhoCoeffs_[coefLabel]
          - molr2*iS.nrhoCoeffs_[coefLabel];
     
        orhoCoeffs_[coefLabel] =
            molr1*orhoCoeffs_[coefLabel]
          - molr2*iS.orhoCoeffs_[coefLabel];

        prhoCoeffs_[coefLabel] =
            molr1*prhoCoeffs_[coefLabel]
          - molr2*iS.prhoCoeffs_[coefLabel];

        qrhoCoeffs_[coefLabel] =
            molr1*qrhoCoeffs_[coefLabel]
          - molr2*iS.qrhoCoeffs_[coefLabel];
 
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::icoSpline<Specie> Foam::operator+
(
    const icoSpline<Specie>& iS1,
    const icoSpline<Specie>& iS2
)
{
    Specie spec = iS1;
     spec += iS2;

    scalar molr1 = iS1.nMoles()/ spec.nMoles();
    scalar molr2 = iS2.nMoles()/ spec.nMoles();

    typename icoSpline<Specie>::coeffArray arhoCoeffs;
    typename icoSpline<Specie>::coeffArray brhoCoeffs;
    typename icoSpline<Specie>::coeffArray crhoCoeffs;
    typename icoSpline<Specie>::coeffArray drhoCoeffs;
    typename icoSpline<Specie>::coeffArray erhoCoeffs;
    typename icoSpline<Specie>::coeffArray frhoCoeffs;
    typename icoSpline<Specie>::coeffArray grhoCoeffs;
    typename icoSpline<Specie>::coeffArray hrhoCoeffs;
    typename icoSpline<Specie>::coeffArray irhoCoeffs;
    typename icoSpline<Specie>::coeffArray jrhoCoeffs;
    typename icoSpline<Specie>::coeffArray krhoCoeffs;
    typename icoSpline<Specie>::coeffArray lrhoCoeffs;
    typename icoSpline<Specie>::coeffArray mrhoCoeffs;
    typename icoSpline<Specie>::coeffArray nrhoCoeffs;
    typename icoSpline<Specie>::coeffArray orhoCoeffs;
    typename icoSpline<Specie>::coeffArray prhoCoeffs;
    typename icoSpline<Specie>::coeffArray qrhoCoeffs;


    for
    (
        register label coefLabel=0;
        coefLabel<icoSpline<Specie>::nCoeffs_;
        coefLabel++
    )
    {
        arhoCoeffs[coefLabel] =
            molr1*iS1.arhoCoeffs_[coefLabel]
          + molr2*iS2.arhoCoeffs_[coefLabel];
        
        brhoCoeffs[coefLabel] =
            molr1*iS1.brhoCoeffs_[coefLabel]
          + molr2*iS2.brhoCoeffs_[coefLabel];

        crhoCoeffs[coefLabel] =
            molr1*iS1.crhoCoeffs_[coefLabel]
          + molr2*iS2.crhoCoeffs_[coefLabel];

        drhoCoeffs[coefLabel] =
            molr1*iS1.drhoCoeffs_[coefLabel]
          + molr2*iS2.drhoCoeffs_[coefLabel];

        erhoCoeffs[coefLabel] =
            molr1*iS1.erhoCoeffs_[coefLabel]
          + molr2*iS2.erhoCoeffs_[coefLabel];

        frhoCoeffs[coefLabel] =
            molr1*iS1.frhoCoeffs_[coefLabel]
          + molr2*iS2.frhoCoeffs_[coefLabel];
        
        grhoCoeffs[coefLabel] =
            molr1*iS1.grhoCoeffs_[coefLabel]
          + molr2*iS2.grhoCoeffs_[coefLabel];

        hrhoCoeffs[coefLabel] =
            molr1*iS1.hrhoCoeffs_[coefLabel]
          + molr2*iS2.hrhoCoeffs_[coefLabel];

        irhoCoeffs[coefLabel] =
            molr1*iS1.irhoCoeffs_[coefLabel]
          + molr2*iS2.irhoCoeffs_[coefLabel];

        jrhoCoeffs[coefLabel] =
            molr1*iS1.jrhoCoeffs_[coefLabel]
          + molr2*iS2.jrhoCoeffs_[coefLabel];

        krhoCoeffs[coefLabel] =
            molr1*iS1.krhoCoeffs_[coefLabel]
          + molr2*iS2.krhoCoeffs_[coefLabel];
        
        lrhoCoeffs[coefLabel] =
            molr1*iS1.lrhoCoeffs_[coefLabel]
          + molr2*iS2.lrhoCoeffs_[coefLabel];

        mrhoCoeffs[coefLabel] =
            molr1*iS1.mrhoCoeffs_[coefLabel]
          + molr2*iS2.mrhoCoeffs_[coefLabel];

        nrhoCoeffs[coefLabel] =
            molr1*iS1.nrhoCoeffs_[coefLabel]
          + molr2*iS2.nrhoCoeffs_[coefLabel];

        orhoCoeffs[coefLabel] =
            molr1*iS1.orhoCoeffs_[coefLabel]
          + molr2*iS2.orhoCoeffs_[coefLabel];

        prhoCoeffs[coefLabel] =
            molr1*iS1.prhoCoeffs_[coefLabel]
          + molr2*iS2.prhoCoeffs_[coefLabel];
        
        qrhoCoeffs[coefLabel] =
            molr1*iS1.qrhoCoeffs_[coefLabel]
          + molr2*iS2.qrhoCoeffs_[coefLabel];

    }

    if
    (
        icoSpline<Specie>::debug
     && notEqual(iS1.T1_, iS2.T1_)
    )
    {
        FatalErrorIn
        (
            "operator+"
            "(const icoSpline<Specie>& iS1,"
            " const icoSpline<Specie>& iS2)"
        )   << "T1 " << iS1.T1_ << " for "
            << (iS1.name().size() ? iS1.name() : "others")
            << " != " << iS2.T1_ << " for "
            << (iS2.name().size() ? iS2.name() : "others")
            << exit(FatalError);
    }

    return icoSpline<Specie>
    (
        spec,
        max(iS1.Tlow_, iS2.Tlow_),
        min(iS1.Thigh_, iS2.Thigh_),
        iS1.T1_,
        iS1.T2_,
        iS1.T3_,
        iS1.T4_,
        iS1.T5_,
        iS1.T6_,
        iS1.T7_,
        iS1.T8_,
        iS1.T9_,
        iS1.T10_,
        iS1.T11_,
        iS1.T12_,
        iS1.T13_,
        iS1.T14_,
        iS1.T15_,
        iS1.T16_,
        arhoCoeffs,
        brhoCoeffs,
        crhoCoeffs,
        drhoCoeffs,
        erhoCoeffs,
        frhoCoeffs,
        grhoCoeffs,
        hrhoCoeffs,
        irhoCoeffs,
        jrhoCoeffs,
        krhoCoeffs,
        lrhoCoeffs,
        mrhoCoeffs,
        nrhoCoeffs,
        orhoCoeffs,
        prhoCoeffs,
        qrhoCoeffs
    );
}


template<class Specie>
inline Foam::icoSpline<Specie> Foam::operator-
(
    const icoSpline<Specie>& iS1,
    const icoSpline<Specie>& iS2
)
{
    Specie spec = iS1;
     spec -= iS2;

    scalar molr1 = iS1.nMoles()/ spec.nMoles();
    scalar molr2 = iS2.nMoles()/ spec.nMoles();

    typename icoSpline<Specie>::coeffArray arhoCoeffs;
    typename icoSpline<Specie>::coeffArray brhoCoeffs;
    typename icoSpline<Specie>::coeffArray crhoCoeffs;
    typename icoSpline<Specie>::coeffArray drhoCoeffs;
    typename icoSpline<Specie>::coeffArray erhoCoeffs;
    typename icoSpline<Specie>::coeffArray frhoCoeffs;
    typename icoSpline<Specie>::coeffArray grhoCoeffs;
    typename icoSpline<Specie>::coeffArray hrhoCoeffs;
    typename icoSpline<Specie>::coeffArray irhoCoeffs;
    typename icoSpline<Specie>::coeffArray jrhoCoeffs;
    typename icoSpline<Specie>::coeffArray krhoCoeffs;
    typename icoSpline<Specie>::coeffArray lrhoCoeffs;
    typename icoSpline<Specie>::coeffArray mrhoCoeffs;
    typename icoSpline<Specie>::coeffArray nrhoCoeffs;
    typename icoSpline<Specie>::coeffArray orhoCoeffs;
    typename icoSpline<Specie>::coeffArray prhoCoeffs;
    typename icoSpline<Specie>::coeffArray qrhoCoeffs;


    for
    (
        register label coefLabel=0;
        coefLabel<icoSpline<Specie>::nCoeffs_;
        coefLabel++
    )
    {
        arhoCoeffs[coefLabel] =
            molr1*iS1.arhoCoeffs_[coefLabel]
          - molr2*iS2.arhoCoeffs_[coefLabel];
        
        brhoCoeffs[coefLabel] =
            molr1*iS1.brhoCoeffs_[coefLabel]
          - molr2*iS2.brhoCoeffs_[coefLabel];

        crhoCoeffs[coefLabel] =
            molr1*iS1.crhoCoeffs_[coefLabel]
          - molr2*iS2.crhoCoeffs_[coefLabel];

        drhoCoeffs[coefLabel] =
            molr1*iS1.drhoCoeffs_[coefLabel]
          - molr2*iS2.drhoCoeffs_[coefLabel];

        erhoCoeffs[coefLabel] =
            molr1*iS1.erhoCoeffs_[coefLabel]
          - molr2*iS2.erhoCoeffs_[coefLabel];

        frhoCoeffs[coefLabel] =
            molr1*iS1.frhoCoeffs_[coefLabel]
          - molr2*iS2.frhoCoeffs_[coefLabel];
        
        grhoCoeffs[coefLabel] =
            molr1*iS1.grhoCoeffs_[coefLabel]
          - molr2*iS2.grhoCoeffs_[coefLabel];

        hrhoCoeffs[coefLabel] =
            molr1*iS1.hrhoCoeffs_[coefLabel]
          - molr2*iS2.hrhoCoeffs_[coefLabel];

        irhoCoeffs[coefLabel] =
            molr1*iS1.irhoCoeffs_[coefLabel]
          - molr2*iS2.irhoCoeffs_[coefLabel];

        jrhoCoeffs[coefLabel] =
            molr1*iS1.jrhoCoeffs_[coefLabel]
          - molr2*iS2.jrhoCoeffs_[coefLabel];

        krhoCoeffs[coefLabel] =
            molr1*iS1.krhoCoeffs_[coefLabel]
          - molr2*iS2.krhoCoeffs_[coefLabel];
        
        lrhoCoeffs[coefLabel] =
            molr1*iS1.lrhoCoeffs_[coefLabel]
          - molr2*iS2.lrhoCoeffs_[coefLabel];

        mrhoCoeffs[coefLabel] =
            molr1*iS1.mrhoCoeffs_[coefLabel]
          - molr2*iS2.mrhoCoeffs_[coefLabel];

        nrhoCoeffs[coefLabel] =
            molr1*iS1.nrhoCoeffs_[coefLabel]
          - molr2*iS2.nrhoCoeffs_[coefLabel];

        orhoCoeffs[coefLabel] =
            molr1*iS1.orhoCoeffs_[coefLabel]
          - molr2*iS2.orhoCoeffs_[coefLabel];

        prhoCoeffs[coefLabel] =
            molr1*iS1.prhoCoeffs_[coefLabel]
          - molr2*iS2.prhoCoeffs_[coefLabel];
        

    }

    if
    (
        icoSpline<Specie>::debug
     && notEqual(iS1.T1_, iS2.T1_)
    )
    {
        FatalErrorIn
        (
            "operator+"
            "(const icoSpline<Specie>& is1,"
            " const icoSpline<Specie>& is2)"
        )   << "T1 " << iS1.T1_ << " for "
            << (iS1.name().size() ? iS1.name() : "others")
            << " != " << iS2.T1_ << " for "
            << (iS2.name().size() ? iS2.name() : "others")
            << exit(FatalError);
    }

    return icoSpline<Specie>
    (
        spec,
        max(iS1.Tlow_, iS2.Tlow_),
        min(iS1.Thigh_, iS2.Thigh_),
        iS1.T1_,
        iS1.T2_,
        iS1.T3_,
        iS1.T4_,
        iS1.T5_,
        iS1.T6_,
        iS1.T7_,
        iS1.T8_,
        iS1.T9_,
        iS1.T10_,
        iS1.T11_,
        iS1.T12_,
        iS1.T13_,
        iS1.T14_,
        iS1.T15_,
        iS1.T16_,
        arhoCoeffs,
        brhoCoeffs,
        crhoCoeffs,
        drhoCoeffs,
        erhoCoeffs,
        frhoCoeffs,
        grhoCoeffs,
        hrhoCoeffs,
        irhoCoeffs,
        jrhoCoeffs,
        krhoCoeffs,
        lrhoCoeffs,
        mrhoCoeffs,
        nrhoCoeffs,
        orhoCoeffs,
        prhoCoeffs,
        qrhoCoeffs
    );
}


template<class Specie>
inline Foam::icoSpline<Specie> Foam::operator*
(
    const scalar s,
    const icoSpline<Specie>& iS
)
{
    return icoSpline<Specie>
    (
        s*static_cast<const Specie&>(iS),
        iS.Tlow_,
        iS.Thigh_,
        iS.T1_,
        iS.T2_,
        iS.T3_,
        iS.T4_,
        iS.T5_,
        iS.T6_,
        iS.T7_,
        iS.T8_,
        iS.T9_,
        iS.T10_,
        iS.T11_,
        iS.T12_,
        iS.T13_,
        iS.T14_,
        iS.T15_,
        iS.T16_,
        iS.aCpCoeffs_,
        iS.bCpCoeffs_,
        iS.cCpCoeffs_,
        iS.dCpCoeffs_,
        iS.eCpCoeffs_,
        iS.fCpCoeffs_,
        iS.gCpCoeffs_,
        iS.hCpCoeffs_,
        iS.iCpCoeffs_,
        iS.jCpCoeffs_,
        iS.kCpCoeffs_,
        iS.lCpCoeffs_,
        iS.mCpCoeffs_,
        iS.nCpCoeffs_,
        iS.oCpCoeffs_,
        iS.pCpCoeffs_,
        iS.qCpCoeffs_
    );
}


template<class Specie>
inline Foam::icoSpline<Specie> Foam::operator==
(
    const icoSpline<Specie>& iS1,
    const icoSpline<Specie>& iS2
)
{
    return iS2 - iS1;
}


// ************************************************************************* //
