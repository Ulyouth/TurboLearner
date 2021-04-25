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

#include "splineTransport.H"
#include "specie.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
inline Foam::splineTransport<Thermo>::splineTransport
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
    const typename splineTransport<Thermo>::coeffArray& amuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& bmuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& cmuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& dmuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& emuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& fmuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& gmuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& hmuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& imuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& jmuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& kmuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& lmuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& mmuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& nmuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& omuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& pmuCoeffs,
    const typename splineTransport<Thermo>::coeffArray& akappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& bkappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& ckappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& dkappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& ekappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& fkappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& gkappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& hkappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& ikappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& jkappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& kkappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& lkappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& mkappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& nkappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& okappaCoeffs,
    const typename splineTransport<Thermo>::coeffArray& pkappaCoeffs
)
:
    Thermo(t),
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
    T15_(T15)

{
    for (register label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        amuCoeffs_[coefLabel] = amuCoeffs[coefLabel];
        bmuCoeffs_[coefLabel] = bmuCoeffs[coefLabel];
        cmuCoeffs_[coefLabel] = cmuCoeffs[coefLabel];
        dmuCoeffs_[coefLabel] = dmuCoeffs[coefLabel];
        emuCoeffs_[coefLabel] = emuCoeffs[coefLabel];
        fmuCoeffs_[coefLabel] = fmuCoeffs[coefLabel];
        gmuCoeffs_[coefLabel] = gmuCoeffs[coefLabel];
        hmuCoeffs_[coefLabel] = hmuCoeffs[coefLabel];
        imuCoeffs_[coefLabel] = imuCoeffs[coefLabel];
        jmuCoeffs_[coefLabel] = jmuCoeffs[coefLabel];
        kmuCoeffs_[coefLabel] = kmuCoeffs[coefLabel];
        lmuCoeffs_[coefLabel] = lmuCoeffs[coefLabel];
        mmuCoeffs_[coefLabel] = mmuCoeffs[coefLabel];
        nmuCoeffs_[coefLabel] = nmuCoeffs[coefLabel];
        omuCoeffs_[coefLabel] = omuCoeffs[coefLabel];
        pmuCoeffs_[coefLabel] = pmuCoeffs[coefLabel];
        akappaCoeffs_[coefLabel] = akappaCoeffs[coefLabel];
        bkappaCoeffs_[coefLabel] = bkappaCoeffs[coefLabel];
        ckappaCoeffs_[coefLabel] = ckappaCoeffs[coefLabel];
        dkappaCoeffs_[coefLabel] = dkappaCoeffs[coefLabel];
        ekappaCoeffs_[coefLabel] = ekappaCoeffs[coefLabel];
        fkappaCoeffs_[coefLabel] = fkappaCoeffs[coefLabel];
        gkappaCoeffs_[coefLabel] = gkappaCoeffs[coefLabel];
        hkappaCoeffs_[coefLabel] = hkappaCoeffs[coefLabel];
        ikappaCoeffs_[coefLabel] = ikappaCoeffs[coefLabel];
        jkappaCoeffs_[coefLabel] = jkappaCoeffs[coefLabel];
        kkappaCoeffs_[coefLabel] = kkappaCoeffs[coefLabel];
        lkappaCoeffs_[coefLabel] = lkappaCoeffs[coefLabel];
        mkappaCoeffs_[coefLabel] = mkappaCoeffs[coefLabel];
        nkappaCoeffs_[coefLabel] = nkappaCoeffs[coefLabel];
        okappaCoeffs_[coefLabel] = okappaCoeffs[coefLabel];
        pkappaCoeffs_[coefLabel] = pkappaCoeffs[coefLabel];
 

    }
}


template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::mucoeffs
(
    const scalar T
) const
{
    if (T < T1_)
    {
        return amuCoeffs_;
    }
    else if (T < T2_)
    {
        return bmuCoeffs_;
    }
    else if (T < T3_)
    {
        return cmuCoeffs_;
    }    
    else if ( T < T4_)
    {
        return dmuCoeffs_;
    }
    else if (T < T5_)
    {
        return emuCoeffs_;
    }
    else if (T < T6_)
    {
        return fmuCoeffs_;
    }
    else if (T < T7_)
    {
        return gmuCoeffs_;
    }    
    else if ( T < T8_)
    {
        return hmuCoeffs_;
    }
    else if (T < T9_)
    {
        return imuCoeffs_;
    }
    else if (T < T10_)
    {
        return jmuCoeffs_;
    }
    else if (T < T11_)
    {
        return kmuCoeffs_;
    }    
    else if ( T < T12_)
    {
        return lmuCoeffs_;
    }
    else if (T < T13_)
    {
        return mmuCoeffs_;
    }
    else if (T < T14_)
    {
        return nmuCoeffs_;
    }
    else if (T < T15_)
    {
        return omuCoeffs_;
    }    
    else
    {
        return pmuCoeffs_;
    }
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::kappacoeffs
(
    const scalar T
) const
{
    if (T < T1_)
    {
        return akappaCoeffs_;
    }
    else if (T < T2_)
    {
        return bkappaCoeffs_;
    }
    else if (T < T3_)
    {
        return ckappaCoeffs_;
    }    
    else if ( T < T4_)
    {
        return dkappaCoeffs_;
    }
    else if (T < T5_)
    {
        return ekappaCoeffs_;
    }
    else if (T < T6_)
    {
        return fkappaCoeffs_;
    }
    else if (T < T7_)
    {
        return gkappaCoeffs_;
    }    
    else if ( T < T8_)
    {
        return hkappaCoeffs_;
    }
    else if (T < T9_)
    {
        return ikappaCoeffs_;
    }
    else if (T < T10_)
    {
        return jkappaCoeffs_;
    }
    else if (T < T11_)
    {
        return kkappaCoeffs_;
    }    
    else if ( T < T12_)
    {
        return lkappaCoeffs_;
    }
    else if (T < T13_)
    {
        return mkappaCoeffs_;
    }
    else if (T < T14_)
    {
        return nkappaCoeffs_;
    }
    else if (T < T15_)
    {
        return okappaCoeffs_;
    }    
    else
    {
        return pkappaCoeffs_;
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::splineTransport<Thermo>::splineTransport
(
    const word& name,
    const splineTransport<Thermo>& spt
)
:
    Thermo(spt),
    Tlow_(spt.Tlow_),
    Thigh_(spt.Thigh_),
    T1_(spt.T1_),
    T2_(spt.T2_),
    T3_(spt.T3_),
    T4_(spt.T4_),
    T5_(spt.T5_),
    T6_(spt.T6_),
    T7_(spt.T7_),
    T8_(spt.T8_),
    T9_(spt.T9_),
    T10_(spt.T10_),
    T11_(spt.T11_),
    T12_(spt.T12_),
    T13_(spt.T13_),
    T14_(spt.T14_),
    T15_(spt.T15_)
{
    for (register label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        amuCoeffs_[coefLabel] = spt.amuCoeffs_[coefLabel];
        bmuCoeffs_[coefLabel] = spt.bmuCoeffs_[coefLabel];
        cmuCoeffs_[coefLabel] = spt.cmuCoeffs_[coefLabel];
        dmuCoeffs_[coefLabel] = spt.dmuCoeffs_[coefLabel];
        emuCoeffs_[coefLabel] = spt.emuCoeffs_[coefLabel];
        fmuCoeffs_[coefLabel] = spt.fmuCoeffs_[coefLabel];
        gmuCoeffs_[coefLabel] = spt.gmuCoeffs_[coefLabel];
        hmuCoeffs_[coefLabel] = spt.hmuCoeffs_[coefLabel];
        imuCoeffs_[coefLabel] = spt.imuCoeffs_[coefLabel];
        jmuCoeffs_[coefLabel] = spt.jmuCoeffs_[coefLabel];
        kmuCoeffs_[coefLabel] = spt.kmuCoeffs_[coefLabel];
        lmuCoeffs_[coefLabel] = spt.lmuCoeffs_[coefLabel];
        mmuCoeffs_[coefLabel] = spt.mmuCoeffs_[coefLabel];
        nmuCoeffs_[coefLabel] = spt.nmuCoeffs_[coefLabel];
        omuCoeffs_[coefLabel] = spt.omuCoeffs_[coefLabel];
        pmuCoeffs_[coefLabel] = spt.pmuCoeffs_[coefLabel];
        akappaCoeffs_[coefLabel] = spt.akappaCoeffs_[coefLabel];
        bkappaCoeffs_[coefLabel] = spt.bkappaCoeffs_[coefLabel];
        ckappaCoeffs_[coefLabel] = spt.ckappaCoeffs_[coefLabel];
        dkappaCoeffs_[coefLabel] = spt.dkappaCoeffs_[coefLabel];
        ekappaCoeffs_[coefLabel] = spt.ekappaCoeffs_[coefLabel];
        fkappaCoeffs_[coefLabel] = spt.fkappaCoeffs_[coefLabel];
        gkappaCoeffs_[coefLabel] = spt.gkappaCoeffs_[coefLabel];
        hkappaCoeffs_[coefLabel] = spt.hkappaCoeffs_[coefLabel];
        ikappaCoeffs_[coefLabel] = spt.ikappaCoeffs_[coefLabel];
        jkappaCoeffs_[coefLabel] = spt.jkappaCoeffs_[coefLabel];
        kkappaCoeffs_[coefLabel] = spt.kkappaCoeffs_[coefLabel];
        lkappaCoeffs_[coefLabel] = spt.lkappaCoeffs_[coefLabel];
        mkappaCoeffs_[coefLabel] = spt.mkappaCoeffs_[coefLabel];
        nkappaCoeffs_[coefLabel] = spt.nkappaCoeffs_[coefLabel];
        okappaCoeffs_[coefLabel] = spt.okappaCoeffs_[coefLabel];
        pkappaCoeffs_[coefLabel] = spt.pkappaCoeffs_[coefLabel];
    }
}

template<class Thermo>
inline Foam::autoPtr<Foam::splineTransport<Thermo> >
Foam::splineTransport<Thermo>::clone() const
{
    return autoPtr<splineTransport<Thermo> >
    (
        new splineTransport<Thermo>(*this)
    );//
}

template<class Thermo>
inline Foam::autoPtr<Foam::splineTransport<Thermo> >
Foam::splineTransport<Thermo>::New(Istream& is)
{
    return autoPtr<splineTransport<Thermo> >
    (
        new splineTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::splineTransport<Thermo> >
Foam::splineTransport<Thermo>::New(const dictionary& dict)
{
    return autoPtr<splineTransport<Thermo> >
    (
        new splineTransport<Thermo>(dict)
    );
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::Thigh() const
{
    return Thigh_;
}


template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T1() const
{
    return T1_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T2() const
{
    return T2_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T3() const
{
    return T3_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T4() const
{
    return T4_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T5() const
{
    return T5_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T6() const
{
    return T6_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T7() const
{
    return T7_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T8() const
{
    return T8_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T9() const
{
    return T9_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T10() const
{
    return T10_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T11() const
{
    return T11_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T12() const
{
    return T12_;
}
template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T13() const
{
    return T13_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T14() const
{
    return T14_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::T15() const
{
    return T15_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::amuCoeffs() const
{
    return amuCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::bmuCoeffs() const
{
    return bmuCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::cmuCoeffs() const
{
    return cmuCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::dmuCoeffs() const
{
    return dmuCoeffs_;
}
template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::emuCoeffs() const
{
    return emuCoeffs_;
}
template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::fmuCoeffs() const
{
    return fmuCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::gmuCoeffs() const
{
    return gmuCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::hmuCoeffs() const
{
    return hmuCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::imuCoeffs() const
{
    return imuCoeffs_;
}
template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::jmuCoeffs() const
{
    return jmuCoeffs_;
}
template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::kmuCoeffs() const
{
    return kmuCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::lmuCoeffs() const
{
    return lmuCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::mmuCoeffs() const
{
    return mmuCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::nmuCoeffs() const
{
    return nmuCoeffs_;
}
template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::omuCoeffs() const
{
    return omuCoeffs_;
}
template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::pmuCoeffs() const
{
    return pmuCoeffs_;
}



template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::akappaCoeffs() const
{
    return akappaCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::bkappaCoeffs() const
{
    return bkappaCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::ckappaCoeffs() const
{
    return ckappaCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::dkappaCoeffs() const
{
    return dkappaCoeffs_;
}
template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::ekappaCoeffs() const
{
    return ekappaCoeffs_;
}
template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::fkappaCoeffs() const
{
    return fkappaCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::gkappaCoeffs() const
{
    return gkappaCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::hkappaCoeffs() const
{
    return hkappaCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::ikappaCoeffs() const
{
    return ikappaCoeffs_;
}
template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::jkappaCoeffs() const
{
    return jkappaCoeffs_;
}
template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::kkappaCoeffs() const
{
    return kkappaCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::lkappaCoeffs() const
{
    return lkappaCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::mkappaCoeffs() const
{
    return mkappaCoeffs_;
}

template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::nkappaCoeffs() const
{
    return nkappaCoeffs_;
}
template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::okappaCoeffs() const
{
    return okappaCoeffs_;
}
template<class Thermo>
inline const typename Foam::splineTransport<Thermo>::coeffArray&
Foam::splineTransport<Thermo>::pkappaCoeffs() const
{
    return pkappaCoeffs_;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::mu
(
    scalar p,
    scalar T
) const
{
     const coeffArray& a = mucoeffs(T);
    return  (a[0]*pow((T-273.15),3)+a[1]*pow((T-273.15),2)+a[2]*pow((T-273.15),1)+a[3])*1e-6;
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::kappa
(
    scalar p,
    scalar T
) const
{
     const coeffArray& a = kappacoeffs(T);
    return  (a[0]*pow((T-273.15),3)+a[1]*pow((T-273.15),2)+a[2]*pow((T-273.15),1)+a[3]);
}

template<class Thermo>
inline Foam::scalar Foam::splineTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = kappacoeffs(T);
    return (a[0]*pow((T-273.15),3)+a[1]*pow((T-273.15),2)+a[2]*pow((T-273.15),1)+a[3])/this->cp(p, T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


template<class Thermo>
inline void Foam::splineTransport<Thermo>::operator+=
(
    const splineTransport<Thermo>& spt
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(spt);

    molr1 /= this->nMoles();
    scalar molr2 = spt.nMoles()/this->nMoles();

    Tlow_ = max(Tlow_, spt.Tlow_);
    Thigh_ = min(Thigh_, spt.Thigh_);

    if (splineTransport<Thermo>::debug && notEqual(T1_, spt.T1_))
    {
        FatalErrorIn
        (
            "splineTransport<Thermo>::operator+="
            "(const splineTransport<Thermo>& spt) const"
        )   << "T1 " << T1_ << " for "
            << (this->name().size() ? this->name() : "others")
            << " != " << spt.T1_ << " for "
            << (spt.name().size() ? spt.name() : "others")
            << exit(FatalError);
    }

    for
    (
        register label coefLabel=0;
        coefLabel<splineTransport<Thermo>::nCoeffs_;
        coefLabel++
    )
    {
        amuCoeffs_[coefLabel] =
            molr1*amuCoeffs_[coefLabel]
          + molr2*spt.amuCoeffs_[coefLabel];

        bmuCoeffs_[coefLabel] =
            molr1*bmuCoeffs_[coefLabel]
          + molr2*spt.bmuCoeffs_[coefLabel];
   
        cmuCoeffs_[coefLabel] =
            molr1*cmuCoeffs_[coefLabel]
          + molr2*spt.cmuCoeffs_[coefLabel];

        dmuCoeffs_[coefLabel] =
            molr1*dmuCoeffs_[coefLabel]
          + molr2*spt.dmuCoeffs_[coefLabel];
     
        emuCoeffs_[coefLabel] =
            molr1*emuCoeffs_[coefLabel]
          + molr2*spt.emuCoeffs_[coefLabel];

        fmuCoeffs_[coefLabel] =
            molr1*fmuCoeffs_[coefLabel]
          + molr2*spt.fmuCoeffs_[coefLabel];

        gmuCoeffs_[coefLabel] =
            molr1*gmuCoeffs_[coefLabel]
          + molr2*spt.gmuCoeffs_[coefLabel];
   
        hmuCoeffs_[coefLabel] =
            molr1*hmuCoeffs_[coefLabel]
          + molr2*spt.hmuCoeffs_[coefLabel];

        imuCoeffs_[coefLabel] =
            molr1*imuCoeffs_[coefLabel]
          + molr2*spt.imuCoeffs_[coefLabel];
     
        jmuCoeffs_[coefLabel] =
            molr1*jmuCoeffs_[coefLabel]
          + molr2*spt.jmuCoeffs_[coefLabel];

        kmuCoeffs_[coefLabel] =
            molr1*kmuCoeffs_[coefLabel]
          + molr2*spt.kmuCoeffs_[coefLabel];

        lmuCoeffs_[coefLabel] =
            molr1*lmuCoeffs_[coefLabel]
          + molr2*spt.lmuCoeffs_[coefLabel];
   
        mmuCoeffs_[coefLabel] =
            molr1*mmuCoeffs_[coefLabel]
          + molr2*spt.mmuCoeffs_[coefLabel];

        nmuCoeffs_[coefLabel] =
            molr1*nmuCoeffs_[coefLabel]
          + molr2*spt.nmuCoeffs_[coefLabel];
     
        omuCoeffs_[coefLabel] =
            molr1*omuCoeffs_[coefLabel]
          + molr2*spt.omuCoeffs_[coefLabel];

        pmuCoeffs_[coefLabel] =
            molr1*pmuCoeffs_[coefLabel]
          + molr2*spt.pmuCoeffs_[coefLabel];



        akappaCoeffs_[coefLabel] =
            molr1*akappaCoeffs_[coefLabel]
          + molr2*spt.akappaCoeffs_[coefLabel];

        bkappaCoeffs_[coefLabel] =
            molr1*bkappaCoeffs_[coefLabel]
          + molr2*spt.bkappaCoeffs_[coefLabel];
   
        ckappaCoeffs_[coefLabel] =
            molr1*ckappaCoeffs_[coefLabel]
          + molr2*spt.ckappaCoeffs_[coefLabel];

        dkappaCoeffs_[coefLabel] =
            molr1*dkappaCoeffs_[coefLabel]
          + molr2*spt.dkappaCoeffs_[coefLabel];
     
        ekappaCoeffs_[coefLabel] =
            molr1*ekappaCoeffs_[coefLabel]
          + molr2*spt.ekappaCoeffs_[coefLabel];

        fkappaCoeffs_[coefLabel] =
            molr1*fkappaCoeffs_[coefLabel]
          + molr2*spt.fkappaCoeffs_[coefLabel];

        gkappaCoeffs_[coefLabel] =
            molr1*gkappaCoeffs_[coefLabel]
          + molr2*spt.gkappaCoeffs_[coefLabel];
   
        hkappaCoeffs_[coefLabel] =
            molr1*hkappaCoeffs_[coefLabel]
          + molr2*spt.hkappaCoeffs_[coefLabel];

        ikappaCoeffs_[coefLabel] =
            molr1*ikappaCoeffs_[coefLabel]
          + molr2*spt.ikappaCoeffs_[coefLabel];
     
        jkappaCoeffs_[coefLabel] =
            molr1*jkappaCoeffs_[coefLabel]
          + molr2*spt.jkappaCoeffs_[coefLabel];

        kkappaCoeffs_[coefLabel] =
            molr1*kkappaCoeffs_[coefLabel]
          + molr2*spt.kkappaCoeffs_[coefLabel];

        lkappaCoeffs_[coefLabel] =
            molr1*lkappaCoeffs_[coefLabel]
          + molr2*spt.lkappaCoeffs_[coefLabel];
   
        mkappaCoeffs_[coefLabel] =
            molr1*mkappaCoeffs_[coefLabel]
          + molr2*spt.mkappaCoeffs_[coefLabel];

        nkappaCoeffs_[coefLabel] =
            molr1*nkappaCoeffs_[coefLabel]
          + molr2*spt.nkappaCoeffs_[coefLabel];
     
        okappaCoeffs_[coefLabel] =
            molr1*okappaCoeffs_[coefLabel]
          + molr2*spt.okappaCoeffs_[coefLabel];

        pkappaCoeffs_[coefLabel] =
            molr1*pkappaCoeffs_[coefLabel]
          + molr2*spt.pkappaCoeffs_[coefLabel];
   
   
    }
}



template<class Thermo>
inline void Foam::splineTransport<Thermo>::operator-=
(
    const splineTransport<Thermo>& spt
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(spt);

    molr1 /= this->nMoles();
    scalar molr2 = spt.nMoles()/this->nMoles();

    Tlow_ = max(Tlow_, spt.Tlow_);
    Thigh_ = min(Thigh_, spt.Thigh_);

    if (splineTransport<Thermo>::debug && notEqual(T1_, spt.T1_))
    {
        FatalErrorIn
        (
            "splineTransport<Thermo>::operator+="
            "(const splineTransport<Thermo>& is) const"
        )   << "T1 " << T1_ << " for "
            << (this->name().size() ? this->name() : "others")
            << " != " << spt.T1_ << " for "
            << (spt.name().size() ? spt.name() : "others")
            << exit(FatalError);
    }

    for
    (
        register label coefLabel=0;
        coefLabel<splineTransport<Thermo>::nCoeffs_;
        coefLabel++
    )
    {
        amuCoeffs_[coefLabel] =
            molr1*amuCoeffs_[coefLabel]
          - molr2*spt.amuCoeffs_[coefLabel];

        bmuCoeffs_[coefLabel] =
            molr1*bmuCoeffs_[coefLabel]
          - molr2*spt.bmuCoeffs_[coefLabel];
   
        cmuCoeffs_[coefLabel] =
            molr1*cmuCoeffs_[coefLabel]
          - molr2*spt.cmuCoeffs_[coefLabel];

        dmuCoeffs_[coefLabel] =
            molr1*dmuCoeffs_[coefLabel]
          - molr2*spt.dmuCoeffs_[coefLabel];
     
        emuCoeffs_[coefLabel] =
            molr1*emuCoeffs_[coefLabel]
          - molr2*spt.emuCoeffs_[coefLabel];

        fmuCoeffs_[coefLabel] =
            molr1*fmuCoeffs_[coefLabel]
          - molr2*spt.fmuCoeffs_[coefLabel];

        gmuCoeffs_[coefLabel] =
            molr1*gmuCoeffs_[coefLabel]
          - molr2*spt.gmuCoeffs_[coefLabel];
   
        hmuCoeffs_[coefLabel] =
            molr1*hmuCoeffs_[coefLabel]
          - molr2*spt.hmuCoeffs_[coefLabel];

        imuCoeffs_[coefLabel] =
            molr1*imuCoeffs_[coefLabel]
          - molr2*spt.imuCoeffs_[coefLabel];
     
        jmuCoeffs_[coefLabel] =
            molr1*jmuCoeffs_[coefLabel]
          - molr2*spt.jmuCoeffs_[coefLabel];

        kmuCoeffs_[coefLabel] =
            molr1*kmuCoeffs_[coefLabel]
          - molr2*spt.kmuCoeffs_[coefLabel];

        lmuCoeffs_[coefLabel] =
            molr1*lmuCoeffs_[coefLabel]
          - molr2*spt.lmuCoeffs_[coefLabel];
   
        mmuCoeffs_[coefLabel] =
            molr1*mmuCoeffs_[coefLabel]
          - molr2*spt.mmuCoeffs_[coefLabel];

        nmuCoeffs_[coefLabel] =
            molr1*nmuCoeffs_[coefLabel]
          - molr2*spt.nmuCoeffs_[coefLabel];
     
        omuCoeffs_[coefLabel] =
            molr1*omuCoeffs_[coefLabel]
          - molr2*spt.omuCoeffs_[coefLabel];

        pmuCoeffs_[coefLabel] =
            molr1*pmuCoeffs_[coefLabel]
          - molr2*spt.pmuCoeffs_[coefLabel];



        akappaCoeffs_[coefLabel] =
            molr1*akappaCoeffs_[coefLabel]
          - molr2*spt.akappaCoeffs_[coefLabel];

        bkappaCoeffs_[coefLabel] =
            molr1*bkappaCoeffs_[coefLabel]
          - molr2*spt.bkappaCoeffs_[coefLabel];
   
        ckappaCoeffs_[coefLabel] =
            molr1*ckappaCoeffs_[coefLabel]
          - molr2*spt.ckappaCoeffs_[coefLabel];

        dkappaCoeffs_[coefLabel] =
            molr1*dkappaCoeffs_[coefLabel]
          - molr2*spt.dkappaCoeffs_[coefLabel];
     
        ekappaCoeffs_[coefLabel] =
            molr1*ekappaCoeffs_[coefLabel]
          - molr2*spt.ekappaCoeffs_[coefLabel];

        fkappaCoeffs_[coefLabel] =
            molr1*fkappaCoeffs_[coefLabel]
          - molr2*spt.fkappaCoeffs_[coefLabel];

        gkappaCoeffs_[coefLabel] =
            molr1*gkappaCoeffs_[coefLabel]
          - molr2*spt.gkappaCoeffs_[coefLabel];
   
        hkappaCoeffs_[coefLabel] =
            molr1*hkappaCoeffs_[coefLabel]
          - molr2*spt.hkappaCoeffs_[coefLabel];

        ikappaCoeffs_[coefLabel] =
            molr1*ikappaCoeffs_[coefLabel]
          - molr2*spt.ikappaCoeffs_[coefLabel];
     
        jkappaCoeffs_[coefLabel] =
            molr1*jkappaCoeffs_[coefLabel]
          - molr2*spt.jkappaCoeffs_[coefLabel];

        kkappaCoeffs_[coefLabel] =
            molr1*kkappaCoeffs_[coefLabel]
          - molr2*spt.kkappaCoeffs_[coefLabel];

        lkappaCoeffs_[coefLabel] =
            molr1*lkappaCoeffs_[coefLabel]
          - molr2*spt.lkappaCoeffs_[coefLabel];
   
        mkappaCoeffs_[coefLabel] =
            molr1*mkappaCoeffs_[coefLabel]
          - molr2*spt.mkappaCoeffs_[coefLabel];

        nkappaCoeffs_[coefLabel] =
            molr1*nkappaCoeffs_[coefLabel]
          - molr2*spt.nkappaCoeffs_[coefLabel];
     
        okappaCoeffs_[coefLabel] =
            molr1*okappaCoeffs_[coefLabel]
          - molr2*spt.okappaCoeffs_[coefLabel];

        pkappaCoeffs_[coefLabel] =
            molr1*pkappaCoeffs_[coefLabel]
          - molr2*spt.pkappaCoeffs_[coefLabel];
 
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::splineTransport<Thermo> Foam::operator+
(
    const splineTransport<Thermo>& spt1,
    const splineTransport<Thermo>& spt2
)
{
    Thermo t = spt1;
           t += spt2;

    scalar molr1 = spt1.nMoles()/ t.nMoles();
    scalar molr2 = spt2.nMoles()/ t.nMoles();

    typename splineTransport<Thermo>::coeffArray amuCoeffs;
    typename splineTransport<Thermo>::coeffArray bmuCoeffs;
    typename splineTransport<Thermo>::coeffArray cmuCoeffs;
    typename splineTransport<Thermo>::coeffArray dmuCoeffs;
    typename splineTransport<Thermo>::coeffArray emuCoeffs;
    typename splineTransport<Thermo>::coeffArray fmuCoeffs;
    typename splineTransport<Thermo>::coeffArray gmuCoeffs;
    typename splineTransport<Thermo>::coeffArray hmuCoeffs;
    typename splineTransport<Thermo>::coeffArray imuCoeffs;
    typename splineTransport<Thermo>::coeffArray jmuCoeffs;
    typename splineTransport<Thermo>::coeffArray kmuCoeffs;
    typename splineTransport<Thermo>::coeffArray lmuCoeffs;
    typename splineTransport<Thermo>::coeffArray mmuCoeffs;
    typename splineTransport<Thermo>::coeffArray nmuCoeffs;
    typename splineTransport<Thermo>::coeffArray omuCoeffs;
    typename splineTransport<Thermo>::coeffArray pmuCoeffs;
    typename splineTransport<Thermo>::coeffArray akappaCoeffs;
    typename splineTransport<Thermo>::coeffArray bkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray ckappaCoeffs;
    typename splineTransport<Thermo>::coeffArray dkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray ekappaCoeffs;
    typename splineTransport<Thermo>::coeffArray fkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray gkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray hkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray ikappaCoeffs;
    typename splineTransport<Thermo>::coeffArray jkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray kkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray lkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray mkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray nkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray okappaCoeffs;
    typename splineTransport<Thermo>::coeffArray pkappaCoeffs;



    for
    (
        register label coefLabel=0;
        coefLabel<splineTransport<Thermo>::nCoeffs_;
        coefLabel++
    )
    {
        amuCoeffs[coefLabel] =
            molr1*spt1.amuCoeffs_[coefLabel]
          + molr2*spt2.amuCoeffs_[coefLabel];

        bmuCoeffs[coefLabel] =
            molr1*spt1.bmuCoeffs_[coefLabel]
          + molr2*spt2.bmuCoeffs_[coefLabel];
   
        cmuCoeffs[coefLabel] =
            molr1*spt1.cmuCoeffs_[coefLabel]
          + molr2*spt2.cmuCoeffs_[coefLabel];

        dmuCoeffs[coefLabel] =
            molr1*spt1.dmuCoeffs_[coefLabel]
          + molr2*spt2.dmuCoeffs_[coefLabel];
     
        emuCoeffs[coefLabel] =
            molr1*spt1.emuCoeffs_[coefLabel]
          + molr2*spt2.emuCoeffs_[coefLabel];

        fmuCoeffs[coefLabel] =
            molr1*spt1.fmuCoeffs_[coefLabel]
          + molr2*spt2.fmuCoeffs_[coefLabel];

        gmuCoeffs[coefLabel] =
            molr1*spt1.gmuCoeffs_[coefLabel]
          + molr2*spt2.gmuCoeffs_[coefLabel];
   
        hmuCoeffs[coefLabel] =
            molr1*spt1.hmuCoeffs_[coefLabel]
          + molr2*spt2.hmuCoeffs_[coefLabel];

        imuCoeffs[coefLabel] =
            molr1*spt1.imuCoeffs_[coefLabel]
          + molr2*spt2.imuCoeffs_[coefLabel];
     
        jmuCoeffs[coefLabel] =
            molr1*spt1.jmuCoeffs_[coefLabel]
          + molr2*spt2.jmuCoeffs_[coefLabel];

        kmuCoeffs[coefLabel] =
            molr1*spt1.kmuCoeffs_[coefLabel]
          + molr2*spt2.kmuCoeffs_[coefLabel];

        lmuCoeffs[coefLabel] =
            molr1*spt1.lmuCoeffs_[coefLabel]
          + molr2*spt2.lmuCoeffs_[coefLabel];
   
        mmuCoeffs[coefLabel] =
            molr1*spt1.mmuCoeffs_[coefLabel]
          + molr2*spt2.mmuCoeffs_[coefLabel];

        nmuCoeffs[coefLabel] =
            molr1*spt1.nmuCoeffs_[coefLabel]
          + molr2*spt2.nmuCoeffs_[coefLabel];
     
        omuCoeffs[coefLabel] =
            molr1*spt1.omuCoeffs_[coefLabel]
          + molr2*spt2.omuCoeffs_[coefLabel];

        pmuCoeffs[coefLabel] =
            molr1*spt1.pmuCoeffs_[coefLabel]
          + molr2*spt2.pmuCoeffs_[coefLabel];



        akappaCoeffs[coefLabel] =
            molr1*spt1.akappaCoeffs_[coefLabel]
          + molr2*spt2.akappaCoeffs_[coefLabel];

        bkappaCoeffs[coefLabel] =
            molr1*spt1.bkappaCoeffs_[coefLabel]
          + molr2*spt2.bkappaCoeffs_[coefLabel];
   
        ckappaCoeffs[coefLabel] =
            molr1*spt1.ckappaCoeffs_[coefLabel]
          + molr2*spt2.ckappaCoeffs_[coefLabel];

        dkappaCoeffs[coefLabel] =
            molr1*spt1.dkappaCoeffs_[coefLabel]
          + molr2*spt2.dkappaCoeffs_[coefLabel];
     
        ekappaCoeffs[coefLabel] =
            molr1*spt1.ekappaCoeffs_[coefLabel]
          + molr2*spt2.ekappaCoeffs_[coefLabel];

        fkappaCoeffs[coefLabel] =
            molr1*spt1.fkappaCoeffs_[coefLabel]
          + molr2*spt2.fkappaCoeffs_[coefLabel];

        gkappaCoeffs[coefLabel] =
            molr1*spt1.gkappaCoeffs_[coefLabel]
          + molr2*spt2.gkappaCoeffs_[coefLabel];
   
        hkappaCoeffs[coefLabel] =
            molr1*spt1.hkappaCoeffs_[coefLabel]
          + molr2*spt2.hkappaCoeffs_[coefLabel];

        ikappaCoeffs[coefLabel] =
            molr1*spt1.ikappaCoeffs_[coefLabel]
          + molr2*spt2.ikappaCoeffs_[coefLabel];
     
        jkappaCoeffs[coefLabel] =
            molr1*spt1.jkappaCoeffs_[coefLabel]
          + molr2*spt2.jkappaCoeffs_[coefLabel];

        kkappaCoeffs[coefLabel] =
            molr1*spt1.kkappaCoeffs_[coefLabel]
          + molr2*spt2.kkappaCoeffs_[coefLabel];

        lkappaCoeffs[coefLabel] =
            molr1*spt1.lkappaCoeffs_[coefLabel]
          + molr2*spt2.lkappaCoeffs_[coefLabel];
   
        mkappaCoeffs[coefLabel] =
            molr1*spt1.mkappaCoeffs_[coefLabel]
          + molr2*spt2.mkappaCoeffs_[coefLabel];

        nkappaCoeffs[coefLabel] =
            molr1*spt1.nkappaCoeffs_[coefLabel]
          + molr2*spt2.nkappaCoeffs_[coefLabel];
     
        okappaCoeffs[coefLabel] =
            molr1*spt1.okappaCoeffs_[coefLabel]
          + molr2*spt2.okappaCoeffs_[coefLabel];

        pkappaCoeffs[coefLabel] =
            molr1*spt1.pkappaCoeffs_[coefLabel]
          + molr2*spt2.pkappaCoeffs_[coefLabel];


    }

    if
    (
        splineTransport<Thermo>::debug
     && notEqual(spt1.T1_, spt2.T1_)
    )
    {
        FatalErrorIn
        (
            "operator+"
            "(const splineTransport<Thermo>& spt1,"
            " const splineTransport<Thermo>& spt2)"
        )   << "T1 " << spt1.T1_ << " for "
            << (spt1.name().size() ? spt1.name() : "others")
            << " != " << spt2.T1_ << " for "
            << (spt2.name().size() ? spt2.name() : "others")
            << exit(FatalError);
    }

    return splineTransport<Thermo>
    (
        t,
        max(spt1.Tlow_, spt2.Tlow_),
        min(spt1.Thigh_, spt2.Thigh_),
        spt1.T1_,
        spt1.T2_,
        spt1.T3_,
        spt1.T4_,
        spt1.T5_,
        spt1.T6_,
        spt1.T7_,
        spt1.T8_,
        spt1.T9_,
        spt1.T10_,
        spt1.T11_,
        spt1.T12_,
        spt1.T13_,
        spt1.T14_,
        spt1.T15_,
        amuCoeffs,
        bmuCoeffs,
        cmuCoeffs,
        dmuCoeffs,
        emuCoeffs,
        fmuCoeffs,
        gmuCoeffs,
        hmuCoeffs,
        imuCoeffs,
        jmuCoeffs,
        kmuCoeffs,
        lmuCoeffs,
        mmuCoeffs,
        nmuCoeffs,
        omuCoeffs,
        pmuCoeffs,
        akappaCoeffs,
        bkappaCoeffs,
        ckappaCoeffs,
        dkappaCoeffs,
        ekappaCoeffs,
        fkappaCoeffs,
        gkappaCoeffs,
        hkappaCoeffs,
        ikappaCoeffs,
        jkappaCoeffs,
        kkappaCoeffs,
        lkappaCoeffs,
        mkappaCoeffs,
        nkappaCoeffs,
        okappaCoeffs,
        pkappaCoeffs

    );
}


template<class Thermo>
inline Foam::splineTransport<Thermo> Foam::operator-
(
    const splineTransport<Thermo>& spt1,
    const splineTransport<Thermo>& spt2
)
{
    Thermo t = spt1;
           t -= spt2;

    scalar molr1 = spt1.nMoles()/ t.nMoles();
    scalar molr2 = spt2.nMoles()/ t.nMoles();

    typename splineTransport<Thermo>::coeffArray amuCoeffs;
    typename splineTransport<Thermo>::coeffArray bmuCoeffs;
    typename splineTransport<Thermo>::coeffArray cmuCoeffs;
    typename splineTransport<Thermo>::coeffArray dmuCoeffs;
    typename splineTransport<Thermo>::coeffArray emuCoeffs;
    typename splineTransport<Thermo>::coeffArray fmuCoeffs;
    typename splineTransport<Thermo>::coeffArray gmuCoeffs;
    typename splineTransport<Thermo>::coeffArray hmuCoeffs;
    typename splineTransport<Thermo>::coeffArray imuCoeffs;
    typename splineTransport<Thermo>::coeffArray jmuCoeffs;
    typename splineTransport<Thermo>::coeffArray kmuCoeffs;
    typename splineTransport<Thermo>::coeffArray lmuCoeffs;
    typename splineTransport<Thermo>::coeffArray mmuCoeffs;
    typename splineTransport<Thermo>::coeffArray nmuCoeffs;
    typename splineTransport<Thermo>::coeffArray omuCoeffs;
    typename splineTransport<Thermo>::coeffArray pmuCoeffs;
    typename splineTransport<Thermo>::coeffArray akappaCoeffs;
    typename splineTransport<Thermo>::coeffArray bkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray ckappaCoeffs;
    typename splineTransport<Thermo>::coeffArray dkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray ekappaCoeffs;
    typename splineTransport<Thermo>::coeffArray fkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray gkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray hkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray ikappaCoeffs;
    typename splineTransport<Thermo>::coeffArray jkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray kkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray lkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray mkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray nkappaCoeffs;
    typename splineTransport<Thermo>::coeffArray okappaCoeffs;
    typename splineTransport<Thermo>::coeffArray pkappaCoeffs;



    for
    (
        register label coefLabel=0;
        coefLabel<splineTransport<Thermo>::nCoeffs_;
        coefLabel++
    )
    {
        amuCoeffs[coefLabel] =
            molr1*spt1.amuCoeffs_[coefLabel]
          - molr2*spt2.amuCoeffs_[coefLabel];

        bmuCoeffs[coefLabel] =
            molr1*spt1.bmuCoeffs_[coefLabel]
          - molr2*spt2.bmuCoeffs_[coefLabel];
   
        cmuCoeffs[coefLabel] =
            molr1*spt1.cmuCoeffs_[coefLabel]
          - molr2*spt2.cmuCoeffs_[coefLabel];

        dmuCoeffs[coefLabel] =
            molr1*spt1.dmuCoeffs_[coefLabel]
          - molr2*spt2.dmuCoeffs_[coefLabel];
     
        emuCoeffs[coefLabel] =
            molr1*spt1.emuCoeffs_[coefLabel]
          - molr2*spt2.emuCoeffs_[coefLabel];

        fmuCoeffs[coefLabel] =
            molr1*spt1.fmuCoeffs_[coefLabel]
          - molr2*spt2.fmuCoeffs_[coefLabel];

        gmuCoeffs[coefLabel] =
            molr1*spt1.gmuCoeffs_[coefLabel]
          - molr2*spt2.gmuCoeffs_[coefLabel];
   
        hmuCoeffs[coefLabel] =
            molr1*spt1.hmuCoeffs_[coefLabel]
          - molr2*spt2.hmuCoeffs_[coefLabel];

        imuCoeffs[coefLabel] =
            molr1*spt1.imuCoeffs_[coefLabel]
          - molr2*spt2.imuCoeffs_[coefLabel];
     
        jmuCoeffs[coefLabel] =
            molr1*spt1.jmuCoeffs_[coefLabel]
          - molr2*spt2.jmuCoeffs_[coefLabel];

        kmuCoeffs[coefLabel] =
            molr1*spt1.kmuCoeffs_[coefLabel]
          - molr2*spt2.kmuCoeffs_[coefLabel];

        lmuCoeffs[coefLabel] =
            molr1*spt1.lmuCoeffs_[coefLabel]
          - molr2*spt2.lmuCoeffs_[coefLabel];
   
        mmuCoeffs[coefLabel] =
            molr1*spt1.mmuCoeffs_[coefLabel]
          - molr2*spt2.mmuCoeffs_[coefLabel];

        nmuCoeffs[coefLabel] =
            molr1*spt1.nmuCoeffs_[coefLabel]
          - molr2*spt2.nmuCoeffs_[coefLabel];
     
        omuCoeffs[coefLabel] =
            molr1*spt1.omuCoeffs_[coefLabel]
          - molr2*spt2.omuCoeffs_[coefLabel];

        pmuCoeffs[coefLabel] =
            molr1*spt1.pmuCoeffs_[coefLabel]
          - molr2*spt2.pmuCoeffs_[coefLabel];



        akappaCoeffs[coefLabel] =
            molr1*spt1.akappaCoeffs_[coefLabel]
          - molr2*spt2.akappaCoeffs_[coefLabel];

        bkappaCoeffs[coefLabel] =
            molr1*spt1.bkappaCoeffs_[coefLabel]
          - molr2*spt2.bkappaCoeffs_[coefLabel];
   
        ckappaCoeffs[coefLabel] =
            molr1*spt1.ckappaCoeffs_[coefLabel]
          - molr2*spt2.ckappaCoeffs_[coefLabel];

        dkappaCoeffs[coefLabel] =
            molr1*spt1.dkappaCoeffs_[coefLabel]
          - molr2*spt2.dkappaCoeffs_[coefLabel];
     
        ekappaCoeffs[coefLabel] =
            molr1*spt1.ekappaCoeffs_[coefLabel]
          - molr2*spt2.ekappaCoeffs_[coefLabel];

        fkappaCoeffs[coefLabel] =
            molr1*spt1.fkappaCoeffs_[coefLabel]
          - molr2*spt2.fkappaCoeffs_[coefLabel];

        gkappaCoeffs[coefLabel] =
            molr1*spt1.gkappaCoeffs_[coefLabel]
          - molr2*spt2.gkappaCoeffs_[coefLabel];
   
        hkappaCoeffs[coefLabel] =
            molr1*spt1.hkappaCoeffs_[coefLabel]
          - molr2*spt2.hkappaCoeffs_[coefLabel];

        ikappaCoeffs[coefLabel] =
            molr1*spt1.ikappaCoeffs_[coefLabel]
          - molr2*spt2.ikappaCoeffs_[coefLabel];
     
        jkappaCoeffs[coefLabel] =
            molr1*spt1.jkappaCoeffs_[coefLabel]
          - molr2*spt2.jkappaCoeffs_[coefLabel];

        kkappaCoeffs[coefLabel] =
            molr1*spt1.kkappaCoeffs_[coefLabel]
          - molr2*spt2.kkappaCoeffs_[coefLabel];

        lkappaCoeffs[coefLabel] =
            molr1*spt1.lkappaCoeffs_[coefLabel]
          - molr2*spt2.lkappaCoeffs_[coefLabel];
   
        mkappaCoeffs[coefLabel] =
            molr1*spt1.mkappaCoeffs_[coefLabel]
          - molr2*spt2.mkappaCoeffs_[coefLabel];

        nkappaCoeffs[coefLabel] =
            molr1*spt1.nkappaCoeffs_[coefLabel]
          - molr2*spt2.nkappaCoeffs_[coefLabel];
     
        okappaCoeffs[coefLabel] =
            molr1*spt1.okappaCoeffs_[coefLabel]
          - molr2*spt2.okappaCoeffs_[coefLabel];

        pkappaCoeffs[coefLabel] =
            molr1*spt1.pkappaCoeffs_[coefLabel]
          - molr2*spt2.pkappaCoeffs_[coefLabel];


    }

    if
    (
        splineTransport<Thermo>::debug
     && notEqual(spt1.T1_, spt2.T1_)
    )
    {
        FatalErrorIn
        (
            "operator+"
            "(const splineTransport<Thermo>& spt1,"
            " const splineTransport<Thermo>& spt2)"
        )   << "T1 " << spt1.T1_ << " for "
            << (spt1.name().size() ? spt1.name() : "others")
            << " != " << spt2.T1_ << " for "
            << (spt2.name().size() ? spt2.name() : "others")
            << exit(FatalError);
    }

    return splineTransport<Thermo>
    (
        t,
        max(spt1.Tlow_, spt2.Tlow_),
        min(spt1.Thigh_, spt2.Thigh_),
        spt1.T1_,
        spt1.T2_,
        spt1.T3_,
        spt1.T4_,
        spt1.T5_,
        spt1.T6_,
        spt1.T7_,
        spt1.T8_,
        spt1.T9_,
        spt1.T10_,
        spt1.T11_,
        spt1.T12_,
        spt1.T13_,
        spt1.T14_,
        spt1.T15_,
        amuCoeffs,
        bmuCoeffs,
        cmuCoeffs,
        dmuCoeffs,
        emuCoeffs,
        fmuCoeffs,
        gmuCoeffs,
        hmuCoeffs,
        imuCoeffs,
        jmuCoeffs,
        kmuCoeffs,
        lmuCoeffs,
        mmuCoeffs,
        nmuCoeffs,
        omuCoeffs,
        pmuCoeffs,
        akappaCoeffs,
        bkappaCoeffs,
        ckappaCoeffs,
        dkappaCoeffs,
        ekappaCoeffs,
        fkappaCoeffs,
        gkappaCoeffs,
        hkappaCoeffs,
        ikappaCoeffs,
        jkappaCoeffs,
        kkappaCoeffs,
        lkappaCoeffs,
        mkappaCoeffs,
        nkappaCoeffs,
        okappaCoeffs,
        pkappaCoeffs

    );
}


template<class Thermo>
inline Foam::splineTransport<Thermo> Foam::operator*
(
    const scalar s,
    const splineTransport<Thermo>& spt
)
{
    return splineTransport<Thermo>
    (
        s*static_cast<const Thermo&>(spt),
        spt.Tlow_,
        spt.Thigh_,
        spt.T1_,
        spt.T2_,
        spt.T3_,
        spt.T4_,
        spt.T5_,
        spt.T6_,
        spt.T7_,
        spt.T8_,
        spt.T9_,
        spt.T10_,
        spt.T11_,
        spt.T12_,
        spt.T13_,
        spt.T14_,
        spt.T15_,

        spt.amuCoeffs_,
        spt.bmuCoeffs_,
        spt.cmuCoeffs_,
        spt.dmuCoeffs_,
        spt.emuCoeffs_,
        spt.fmuCoeffs_,
        spt.gmuCoeffs_,
        spt.hmuCoeffs_,
        spt.imuCoeffs_,
        spt.jmuCoeffs_,
        spt.kmuCoeffs_,
        spt.lmuCoeffs_,
        spt.mmuCoeffs_,
        spt.nmuCoeffs_,
        spt.omuCoeffs_,
        spt.pmuCoeffs_,
        spt.akappaCoeffs_,
        spt.bkappaCoeffs_,
        spt.ckappaCoeffs_,
        spt.dkappaCoeffs_,
        spt.ekappaCoeffs_,
        spt.fkappaCoeffs_,
        spt.gkappaCoeffs_,
        spt.hkappaCoeffs_,
        spt.ikappaCoeffs_,
        spt.jkappaCoeffs_,
        spt.kkappaCoeffs_,
        spt.lkappaCoeffs_,
        spt.mkappaCoeffs_,
        spt.nkappaCoeffs_,
        spt.okappaCoeffs_,
        spt.pkappaCoeffs_

    );
}


template<class Thermo>
inline Foam::splineTransport<Thermo> Foam::operator==
(
    const splineTransport<Thermo>& spt1,
    const splineTransport<Thermo>& spt2
)
{
    return spt2 - spt1;
}


// ************************************************************************* //
