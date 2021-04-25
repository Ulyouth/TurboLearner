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

#include "splineThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::splineThermo<EquationOfState>::splineThermo
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
    const typename splineThermo<EquationOfState>::coeffArray& a0CpCoeffs,   
    const typename splineThermo<EquationOfState>::coeffArray& aCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& bCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& cCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& dCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& eCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& fCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& gCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& hCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& iCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& jCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& kCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& lCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& mCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& nCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& oCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& pCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& qCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& rCpCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& a0HsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& aHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& bHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& cHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& dHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& eHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& fHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& gHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& hHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& iHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& jHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& kHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& lHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& mHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& nHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& oHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& pHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& qHsCoeffs,
    const typename splineThermo<EquationOfState>::coeffArray& rHsCoeffs
)
:
    EquationOfState(st),
    Tlow_(Tlow),
    Thigh_(Thigh),
    T0_(T0),
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
    T16_(T16),
    T17_(T17)

{
    for (register label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        a0CpCoeffs_[coefLabel] = a0CpCoeffs[coefLabel];
        aCpCoeffs_[coefLabel] = aCpCoeffs[coefLabel];
        bCpCoeffs_[coefLabel] = bCpCoeffs[coefLabel];
        cCpCoeffs_[coefLabel] = cCpCoeffs[coefLabel];
        dCpCoeffs_[coefLabel] = dCpCoeffs[coefLabel];
        eCpCoeffs_[coefLabel] = eCpCoeffs[coefLabel];
        fCpCoeffs_[coefLabel] = fCpCoeffs[coefLabel];
        gCpCoeffs_[coefLabel] = gCpCoeffs[coefLabel];
        hCpCoeffs_[coefLabel] = hCpCoeffs[coefLabel];
        iCpCoeffs_[coefLabel] = iCpCoeffs[coefLabel];
        jCpCoeffs_[coefLabel] = jCpCoeffs[coefLabel];
        kCpCoeffs_[coefLabel] = kCpCoeffs[coefLabel];
        lCpCoeffs_[coefLabel] = lCpCoeffs[coefLabel];
        mCpCoeffs_[coefLabel] = mCpCoeffs[coefLabel];
        nCpCoeffs_[coefLabel] = nCpCoeffs[coefLabel];
        oCpCoeffs_[coefLabel] = oCpCoeffs[coefLabel];
        pCpCoeffs_[coefLabel] = pCpCoeffs[coefLabel];
        qCpCoeffs_[coefLabel] = qCpCoeffs[coefLabel];
        rCpCoeffs_[coefLabel] = rCpCoeffs[coefLabel];
        a0HsCoeffs_[coefLabel] = a0HsCoeffs[coefLabel];
        aHsCoeffs_[coefLabel] = aHsCoeffs[coefLabel];
        bHsCoeffs_[coefLabel] = bHsCoeffs[coefLabel];
        cHsCoeffs_[coefLabel] = cHsCoeffs[coefLabel];
        dHsCoeffs_[coefLabel] = dHsCoeffs[coefLabel];
        eHsCoeffs_[coefLabel] = eHsCoeffs[coefLabel];
        fHsCoeffs_[coefLabel] = fHsCoeffs[coefLabel];
        gHsCoeffs_[coefLabel] = gHsCoeffs[coefLabel];
        hHsCoeffs_[coefLabel] = hHsCoeffs[coefLabel];
        iHsCoeffs_[coefLabel] = iHsCoeffs[coefLabel];
        jHsCoeffs_[coefLabel] = jHsCoeffs[coefLabel];
        kHsCoeffs_[coefLabel] = kHsCoeffs[coefLabel];
        lHsCoeffs_[coefLabel] = lHsCoeffs[coefLabel];
        mHsCoeffs_[coefLabel] = mHsCoeffs[coefLabel];
        nHsCoeffs_[coefLabel] = nHsCoeffs[coefLabel];
        oHsCoeffs_[coefLabel] = oHsCoeffs[coefLabel];
        pHsCoeffs_[coefLabel] = pHsCoeffs[coefLabel];
        qHsCoeffs_[coefLabel] = qHsCoeffs[coefLabel];
        rHsCoeffs_[coefLabel] = rHsCoeffs[coefLabel];
        
    }
}


template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::Cpcoeffs
(
    const scalar T
) const
{
    if (T < T0_)
    {
        return a0CpCoeffs_;
    }

    else if (T < T1_)
    {
        return aCpCoeffs_;
    }
    else if (T < T2_)
    {
        return bCpCoeffs_;
    }
    else if (T < T3_)
    {
        return cCpCoeffs_;
    }    
    else if ( T < T4_)
    {
        return dCpCoeffs_;
    }
    else if (T < T5_)
    {
        return eCpCoeffs_;
    }
    else if (T < T6_)
    {
        return fCpCoeffs_;
    }    
    else if ( T < T7_)
    {
        return gCpCoeffs_;
    }
    else if (T < T8_)
    {
        return hCpCoeffs_;
    }
    else if (T < T9_)
    {
        return iCpCoeffs_;
    }    
    else if ( T < T10_)
    {
        return jCpCoeffs_;
    }
    else if (T < T11_)
    {
        return kCpCoeffs_;
    }
    else if (T < T12_)
    {
        return lCpCoeffs_;
    }    
    else if ( T < T13_)
    {
        return mCpCoeffs_;
    }
    else if ( T < T14_)
    {
        return nCpCoeffs_;
    }
    else if (T < T15_)
    {
        return oCpCoeffs_;
    }
    else if (T < T16_)
    {
        return pCpCoeffs_;
    }

    else if (T < T17_)
    {
        return qCpCoeffs_;
    }      
    else
    {
        return rCpCoeffs_;
    }
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::Hscoeffs
(
    const scalar T
) const
{
    if (T < T0_)
    {
        return a0HsCoeffs_;
    }
    else if (T < T1_)
    {
        return aHsCoeffs_;
    }
    else if (T < T2_)
    {
        return bHsCoeffs_;
    }
    else if (T < T3_)
    {
        return cHsCoeffs_;
    }    
    else if ( T < T4_)
    {
        return dHsCoeffs_;
    }
    else if (T < T5_)
    {
        return eHsCoeffs_;
    }
    else if (T < T6_)
    {
        return fHsCoeffs_;
    }    
    else if ( T < T7_)
    {
        return gHsCoeffs_;
    }
    else if (T < T8_)
    {
        return hHsCoeffs_;
    }
    else if (T < T9_)
    {
        return iHsCoeffs_;
    }    
    else if ( T < T10_)
    {
        return jHsCoeffs_;
    }
    else if (T < T11_)
    {
        return kHsCoeffs_;
    }
    else if (T < T12_)
    {
        return lHsCoeffs_;
    }    
    else if ( T < T13_)
    {
        return mHsCoeffs_;
    }
    else if ( T < T14_)
    {
        return nHsCoeffs_;
    }
    else if (T < T15_)
    {
        return oHsCoeffs_;
    }
    else if (T < T16_)
    {
        return pHsCoeffs_;
    }  
    else if (T < T17_)
    {
        return qHsCoeffs_;
    }    
    else
    {
        return rHsCoeffs_;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::splineThermo<EquationOfState>::splineThermo
(
    const word& name,
    const splineThermo& sth
)
:
    EquationOfState(name, sth),
    Tlow_(sth.Tlow_),
    Thigh_(sth.Thigh_),
    T0_(sth.T0_),
    T1_(sth.T1_),
    T2_(sth.T2_),
    T3_(sth.T3_),
    T4_(sth.T4_),
    T5_(sth.T5_),
    T6_(sth.T6_),
    T7_(sth.T7_),
    T8_(sth.T8_),
    T9_(sth.T9_),
    T10_(sth.T10_),
    T11_(sth.T11_),
    T12_(sth.T12_),
    T13_(sth.T13_),
    T14_(sth.T14_),
    T15_(sth.T15_),
    T16_(sth.T16_),
    T17_(sth.T17_)
{
    for (register label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        a0CpCoeffs_[coefLabel] = sth.a0CpCoeffs_[coefLabel];
        aCpCoeffs_[coefLabel] = sth.aCpCoeffs_[coefLabel];
        bCpCoeffs_[coefLabel] = sth.bCpCoeffs_[coefLabel];
        cCpCoeffs_[coefLabel] = sth.cCpCoeffs_[coefLabel];
        dCpCoeffs_[coefLabel] = sth.dCpCoeffs_[coefLabel];
        eCpCoeffs_[coefLabel] = sth.eCpCoeffs_[coefLabel];
        fCpCoeffs_[coefLabel] = sth.fCpCoeffs_[coefLabel];
        gCpCoeffs_[coefLabel] = sth.gCpCoeffs_[coefLabel];
        hCpCoeffs_[coefLabel] = sth.hCpCoeffs_[coefLabel];
        iCpCoeffs_[coefLabel] = sth.iCpCoeffs_[coefLabel];
        jCpCoeffs_[coefLabel] = sth.jCpCoeffs_[coefLabel];
        kCpCoeffs_[coefLabel] = sth.kCpCoeffs_[coefLabel];
        lCpCoeffs_[coefLabel] = sth.lCpCoeffs_[coefLabel];
        mCpCoeffs_[coefLabel] = sth.mCpCoeffs_[coefLabel];
        nCpCoeffs_[coefLabel] = sth.nCpCoeffs_[coefLabel];
        oCpCoeffs_[coefLabel] = sth.oCpCoeffs_[coefLabel];
        pCpCoeffs_[coefLabel] = sth.pCpCoeffs_[coefLabel];
        qCpCoeffs_[coefLabel] = sth.qCpCoeffs_[coefLabel];
        rCpCoeffs_[coefLabel] = sth.rCpCoeffs_[coefLabel];
        a0HsCoeffs_[coefLabel] = sth.a0HsCoeffs_[coefLabel];
        aHsCoeffs_[coefLabel] = sth.aHsCoeffs_[coefLabel];
        bHsCoeffs_[coefLabel] = sth.bHsCoeffs_[coefLabel];
        cHsCoeffs_[coefLabel] = sth.cHsCoeffs_[coefLabel];
        dHsCoeffs_[coefLabel] = sth.dHsCoeffs_[coefLabel];
        eHsCoeffs_[coefLabel] = sth.eHsCoeffs_[coefLabel];
        fHsCoeffs_[coefLabel] = sth.fHsCoeffs_[coefLabel];
        gHsCoeffs_[coefLabel] = sth.gHsCoeffs_[coefLabel];
        hHsCoeffs_[coefLabel] = sth.hHsCoeffs_[coefLabel];
        iHsCoeffs_[coefLabel] = sth.iHsCoeffs_[coefLabel];
        jHsCoeffs_[coefLabel] = sth.jHsCoeffs_[coefLabel];
        kHsCoeffs_[coefLabel] = sth.kHsCoeffs_[coefLabel];
        lHsCoeffs_[coefLabel] = sth.lHsCoeffs_[coefLabel];
        mHsCoeffs_[coefLabel] = sth.mHsCoeffs_[coefLabel];
        nHsCoeffs_[coefLabel] = sth.nHsCoeffs_[coefLabel];
        oHsCoeffs_[coefLabel] = sth.oHsCoeffs_[coefLabel];
        pHsCoeffs_[coefLabel] = sth.pHsCoeffs_[coefLabel];
        qHsCoeffs_[coefLabel] = sth.qHsCoeffs_[coefLabel];
        rHsCoeffs_[coefLabel] = sth.rHsCoeffs_[coefLabel];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    //if (T < Tlow_ || T > Thigh_)
    //{
    //    WarningIn
    //    (
    //        "splineThermo<EquationOfState>::limit(const scalar T) const"
    //    )   << "attempt to use splineThermo<EquationOfState>"
    //           " out of temperature range "
    //        << Tlow_ << " -> " << Thigh_ << ";  T = " << T
    //        << endl;
    //
    //    return min(max(T, Tlow_), Thigh_);
    //}
    //else
    //{
        return T;
    //}
}


template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::Tlow() const
{
    return Tlow_;
}


template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::Thigh() const
{
    return Thigh_;
}


template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T0() const
{
    return T0_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T1() const
{
    return T1_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T2() const
{
    return T2_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T3() const
{
    return T3_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T4() const
{
    return T4_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T5() const
{
    return T5_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T6() const
{
    return T6_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T7() const
{
    return T7_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T8() const
{
    return T8_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T9() const
{
    return T9_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T10() const
{
    return T10_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T11() const
{
    return T11_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T12() const
{
    return T12_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T13() const
{
    return T13_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T14() const
{
    return T14_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T15() const
{
    return T15_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T16() const
{
    return T16_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::T17() const
{
    return T17_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::a0CpCoeffs() const
{
    return a0CpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::aCpCoeffs() const
{
    return aCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::bCpCoeffs() const
{
    return bCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::cCpCoeffs() const
{
    return cCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::dCpCoeffs() const
{
    return dCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::eCpCoeffs() const
{
    return eCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::fCpCoeffs() const
{
    return fCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::gCpCoeffs() const
{
    return gCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::hCpCoeffs() const
{
    return hCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::iCpCoeffs() const
{
    return iCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::jCpCoeffs() const
{
    return jCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::kCpCoeffs() const
{
    return kCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::lCpCoeffs() const
{
    return lCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::mCpCoeffs() const
{
    return mCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::nCpCoeffs() const
{
    return nCpCoeffs_;
}
template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::oCpCoeffs() const
{
    return oCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::pCpCoeffs() const
{
    return pCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::qCpCoeffs() const
{
    return qCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::rCpCoeffs() const
{
    return rCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::a0HsCoeffs() const
{
    return a0HsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::aHsCoeffs() const
{
    return aHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::bHsCoeffs() const
{
    return bHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::cHsCoeffs() const
{
    return cHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::dHsCoeffs() const
{
    return dHsCoeffs_;
}
template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::eHsCoeffs() const
{
    return eHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::fHsCoeffs() const
{
    return fHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::gHsCoeffs() const
{
    return gHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::hHsCoeffs() const
{
    return hHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::iHsCoeffs() const
{
    return iHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::jHsCoeffs() const
{
    return jHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::kHsCoeffs() const
{
    return kHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::lHsCoeffs() const
{
    return lHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::mHsCoeffs() const
{
    return mHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::nHsCoeffs() const
{
    return nHsCoeffs_;
}
template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::oHsCoeffs() const
{
    return oHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::pHsCoeffs() const
{
    return pHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::qHsCoeffs() const
{
    return qHsCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::splineThermo<EquationOfState>::coeffArray&
Foam::splineThermo<EquationOfState>::rHsCoeffs() const
{
    return rHsCoeffs_;
}

template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = Cpcoeffs(T);
    return (a[0]*pow((T-273.15),3)+a[1]*pow((T-273.15),2)+a[2]*pow((T-273.15),1)+a[3]); 
}


template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::Ha
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = Hscoeffs(T);
    return (a[0]*pow((T-273.15),3)+a[1]*pow((T-273.15),2)+a[2]*pow((T-273.15),1)+a[3])*1000;
}


template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::Hs
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = Hscoeffs(T);
    return Ha(p, T) - Hc();//return (a[0]*pow((T-273.15),3)+a[1]*pow((T-273.15),2)+a[2]*pow((T-273.15),1)+a[3])*1000;//
}


template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::Hc() const

{

    return 0;
}



template<class EquationOfState>
inline Foam::scalar Foam::splineThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
        notImplemented
    (
        "scalar hConstThermo<equationOfState>::"
        "s(const scalar p, const scalar T) const"
    );
    return T;
    //const coeffArray& a = coeffs(T);
    //return
    //this->RR*
    //(
    //    (((a[4]/4.0*T + a[3]/3.0)*T + a[2]/2.0)*T + a[1])*T + a[0]*::log(T)
    //  + a[6]
    //);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


template<class EquationOfState>
inline void Foam::splineThermo<EquationOfState>::operator+=
(
    const splineThermo<EquationOfState>& sth
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator+=(sth);

    molr1 /= this->nMoles();
    scalar molr2 = sth.nMoles()/this->nMoles();

    Tlow_ = max(Tlow_, sth.Tlow_);
    Thigh_ = min(Thigh_, sth.Thigh_);

    if (splineThermo<EquationOfState>::debug && notEqual(T0_, sth.T0_))
    {
        FatalErrorIn
        (
            "splineThermo<EquationOfState>::operator+="
            "(const splineThermo<EquationOfState>& st) const"
        )   << "T0 " << T0_ << " for "
            << (this->name().size() ? this->name() : "others")
            << " != " << sth.T0_ << " for "
            << (sth.name().size() ? sth.name() : "others")
            << exit(FatalError);
    }

    for
    (
        register label coefLabel=0;
        coefLabel<splineThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {

        a0CpCoeffs_[coefLabel] =
            molr1*a0CpCoeffs_[coefLabel]
          + molr2*sth.a0CpCoeffs_[coefLabel];

        aCpCoeffs_[coefLabel] =
            molr1*aCpCoeffs_[coefLabel]
          + molr2*sth.aCpCoeffs_[coefLabel];

        bCpCoeffs_[coefLabel] =
            molr1*bCpCoeffs_[coefLabel]
          + molr2*sth.bCpCoeffs_[coefLabel];
   
        cCpCoeffs_[coefLabel] =
            molr1*cCpCoeffs_[coefLabel]
          + molr2*sth.cCpCoeffs_[coefLabel];

        dCpCoeffs_[coefLabel] =
            molr1*dCpCoeffs_[coefLabel]
          + molr2*sth.dCpCoeffs_[coefLabel];
     
        eCpCoeffs_[coefLabel] =
            molr1*eCpCoeffs_[coefLabel]
          + molr2*sth.eCpCoeffs_[coefLabel];

        fCpCoeffs_[coefLabel] =
            molr1*fCpCoeffs_[coefLabel]
          + molr2*sth.fCpCoeffs_[coefLabel];

        gCpCoeffs_[coefLabel] =
            molr1*gCpCoeffs_[coefLabel]
          + molr2*sth.gCpCoeffs_[coefLabel];
   
        hCpCoeffs_[coefLabel] =
            molr1*hCpCoeffs_[coefLabel]
          + molr2*sth.hCpCoeffs_[coefLabel];

        iCpCoeffs_[coefLabel] =
            molr1*iCpCoeffs_[coefLabel]
          + molr2*sth.iCpCoeffs_[coefLabel];
     
        jCpCoeffs_[coefLabel] =
            molr1*jCpCoeffs_[coefLabel]
          + molr2*sth.jCpCoeffs_[coefLabel];

        kCpCoeffs_[coefLabel] =
            molr1*kCpCoeffs_[coefLabel]
          + molr2*sth.kCpCoeffs_[coefLabel];

        lCpCoeffs_[coefLabel] =
            molr1*lCpCoeffs_[coefLabel]
          + molr2*sth.lCpCoeffs_[coefLabel];
   
        mCpCoeffs_[coefLabel] =
            molr1*mCpCoeffs_[coefLabel]
          + molr2*sth.mCpCoeffs_[coefLabel];

        nCpCoeffs_[coefLabel] =
            molr1*nCpCoeffs_[coefLabel]
          + molr2*sth.nCpCoeffs_[coefLabel];
     
        oCpCoeffs_[coefLabel] =
            molr1*oCpCoeffs_[coefLabel]
          + molr2*sth.oCpCoeffs_[coefLabel];

        pCpCoeffs_[coefLabel] =
            molr1*pCpCoeffs_[coefLabel]
          + molr2*sth.pCpCoeffs_[coefLabel];
     
        qCpCoeffs_[coefLabel] =
            molr1*qCpCoeffs_[coefLabel]
          + molr2*sth.qCpCoeffs_[coefLabel];

        rCpCoeffs_[coefLabel] =
            molr1*rCpCoeffs_[coefLabel]
          + molr2*sth.rCpCoeffs_[coefLabel];

        a0HsCoeffs_[coefLabel] =
            molr1*a0HsCoeffs_[coefLabel]
          + molr2*sth.a0HsCoeffs_[coefLabel];

        aHsCoeffs_[coefLabel] =
            molr1*aHsCoeffs_[coefLabel]
          + molr2*sth.aHsCoeffs_[coefLabel];

        bHsCoeffs_[coefLabel] =
            molr1*bHsCoeffs_[coefLabel]
          + molr2*sth.bHsCoeffs_[coefLabel];
   
        cHsCoeffs_[coefLabel] =
            molr1*cHsCoeffs_[coefLabel]
          + molr2*sth.cHsCoeffs_[coefLabel];

        dHsCoeffs_[coefLabel] =
            molr1*dHsCoeffs_[coefLabel]
          + molr2*sth.dHsCoeffs_[coefLabel];
     
        eHsCoeffs_[coefLabel] =
            molr1*eHsCoeffs_[coefLabel]
          + molr2*sth.eHsCoeffs_[coefLabel];

        fHsCoeffs_[coefLabel] =
            molr1*fHsCoeffs_[coefLabel]
          + molr2*sth.fHsCoeffs_[coefLabel];

        gHsCoeffs_[coefLabel] =
            molr1*gHsCoeffs_[coefLabel]
          + molr2*sth.gHsCoeffs_[coefLabel];
   
        hHsCoeffs_[coefLabel] =
            molr1*hHsCoeffs_[coefLabel]
          + molr2*sth.hHsCoeffs_[coefLabel];

        iHsCoeffs_[coefLabel] =
            molr1*iHsCoeffs_[coefLabel]
          + molr2*sth.iHsCoeffs_[coefLabel];
     
        jHsCoeffs_[coefLabel] =
            molr1*jHsCoeffs_[coefLabel]
          + molr2*sth.jHsCoeffs_[coefLabel];

        kHsCoeffs_[coefLabel] =
            molr1*kHsCoeffs_[coefLabel]
          + molr2*sth.kHsCoeffs_[coefLabel];

        lHsCoeffs_[coefLabel] =
            molr1*lHsCoeffs_[coefLabel]
          + molr2*sth.lHsCoeffs_[coefLabel];
   
        mHsCoeffs_[coefLabel] =
            molr1*mHsCoeffs_[coefLabel]
          + molr2*sth.mHsCoeffs_[coefLabel];

        nHsCoeffs_[coefLabel] =
            molr1*nHsCoeffs_[coefLabel]
          + molr2*sth.nHsCoeffs_[coefLabel];
     
        oHsCoeffs_[coefLabel] =
            molr1*oHsCoeffs_[coefLabel]
          + molr2*sth.oHsCoeffs_[coefLabel];

        pHsCoeffs_[coefLabel] =
            molr1*pHsCoeffs_[coefLabel]
          + molr2*sth.pHsCoeffs_[coefLabel];
     
        qHsCoeffs_[coefLabel] =
            molr1*qHsCoeffs_[coefLabel]
          + molr2*sth.qHsCoeffs_[coefLabel];

        rHsCoeffs_[coefLabel] =
            molr1*rHsCoeffs_[coefLabel]
          + molr2*sth.rHsCoeffs_[coefLabel];
    }
}


template<class EquationOfState>
inline void Foam::splineThermo<EquationOfState>::operator-=
(
    const splineThermo<EquationOfState>& sth
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator-=(sth);

    molr1 /= this->nMoles();
    scalar molr2 = sth.nMoles()/this->nMoles();

    Tlow_ = max(Tlow_, sth.Tlow_);
    Thigh_ = min(Thigh_, sth.Thigh_);

    if (splineThermo<EquationOfState>::debug && notEqual(T0_, sth.T0_))
    {
        FatalErrorIn
        (
            "splineThermo<EquationOfState>::operator-="
            "(const splineThermo<EquationOfState>& st) const"
        )   << "T0 " << T0_ << " for "
            << (this->name().size() ? this->name() : "others")
            << " != " << sth.T0_ << " for "
            << (sth.name().size() ? sth.name() : "others")
            << exit(FatalError);
    }

    for
    (
        register label coefLabel=0;
        coefLabel<splineThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {

        a0CpCoeffs_[coefLabel] =
            molr1*a0CpCoeffs_[coefLabel]
          - molr2*sth.a0CpCoeffs_[coefLabel];

        aCpCoeffs_[coefLabel] =
            molr1*aCpCoeffs_[coefLabel]
          - molr2*sth.aCpCoeffs_[coefLabel];

        bCpCoeffs_[coefLabel] =
            molr1*bCpCoeffs_[coefLabel]
          - molr2*sth.bCpCoeffs_[coefLabel];

        cCpCoeffs_[coefLabel] =
            molr1*cCpCoeffs_[coefLabel]
          - molr2*sth.cCpCoeffs_[coefLabel];

        dCpCoeffs_[coefLabel] =
            molr1*dCpCoeffs_[coefLabel]
          - molr2*sth.dCpCoeffs_[coefLabel];

        eCpCoeffs_[coefLabel] =
            molr1*eCpCoeffs_[coefLabel]
          - molr2*sth.eCpCoeffs_[coefLabel];

        fCpCoeffs_[coefLabel] =
            molr1*fCpCoeffs_[coefLabel]
          - molr2*sth.fCpCoeffs_[coefLabel];

        gCpCoeffs_[coefLabel] =
            molr1*gCpCoeffs_[coefLabel]
          - molr2*sth.gCpCoeffs_[coefLabel];

        hCpCoeffs_[coefLabel] =
            molr1*hCpCoeffs_[coefLabel]
          - molr2*sth.hCpCoeffs_[coefLabel];

        iCpCoeffs_[coefLabel] =
            molr1*iCpCoeffs_[coefLabel]
          - molr2*sth.iCpCoeffs_[coefLabel];

        jCpCoeffs_[coefLabel] =
            molr1*jCpCoeffs_[coefLabel]
          - molr2*sth.jCpCoeffs_[coefLabel];

        kCpCoeffs_[coefLabel] =
            molr1*kCpCoeffs_[coefLabel]
          - molr2*sth.kCpCoeffs_[coefLabel];

        lCpCoeffs_[coefLabel] =
            molr1*lCpCoeffs_[coefLabel]
          - molr2*sth.lCpCoeffs_[coefLabel];

        mCpCoeffs_[coefLabel] =
            molr1*mCpCoeffs_[coefLabel]
          - molr2*sth.mCpCoeffs_[coefLabel];

        nCpCoeffs_[coefLabel] =
            molr1*nCpCoeffs_[coefLabel]
          - molr2*sth.nCpCoeffs_[coefLabel];

        oCpCoeffs_[coefLabel] =
            molr1*oCpCoeffs_[coefLabel]
          - molr2*sth.oCpCoeffs_[coefLabel];

        pCpCoeffs_[coefLabel] =
            molr1*pCpCoeffs_[coefLabel]
          - molr2*sth.pCpCoeffs_[coefLabel];

        qCpCoeffs_[coefLabel] =
            molr1*qCpCoeffs_[coefLabel]
          - molr2*sth.qCpCoeffs_[coefLabel];

        rCpCoeffs_[coefLabel] =
            molr1*rCpCoeffs_[coefLabel]
          - molr2*sth.rCpCoeffs_[coefLabel];

        a0HsCoeffs_[coefLabel] =
            molr1*a0HsCoeffs_[coefLabel]
          - molr2*sth.a0HsCoeffs_[coefLabel];

        aHsCoeffs_[coefLabel] =
            molr1*aHsCoeffs_[coefLabel]
          - molr2*sth.aHsCoeffs_[coefLabel];

        bHsCoeffs_[coefLabel] =
            molr1*bHsCoeffs_[coefLabel]
          - molr2*sth.bHsCoeffs_[coefLabel];
   
        cHsCoeffs_[coefLabel] =
            molr1*cHsCoeffs_[coefLabel]
          - molr2*sth.cHsCoeffs_[coefLabel];

        dHsCoeffs_[coefLabel] =
            molr1*dHsCoeffs_[coefLabel]
          - molr2*sth.dHsCoeffs_[coefLabel];
     
        eHsCoeffs_[coefLabel] =
            molr1*eHsCoeffs_[coefLabel]
          - molr2*sth.eHsCoeffs_[coefLabel];

        fHsCoeffs_[coefLabel] =
            molr1*fHsCoeffs_[coefLabel]
          - molr2*sth.fHsCoeffs_[coefLabel];

        gHsCoeffs_[coefLabel] =
            molr1*gHsCoeffs_[coefLabel]
          - molr2*sth.gHsCoeffs_[coefLabel];
   
        hHsCoeffs_[coefLabel] =
            molr1*hHsCoeffs_[coefLabel]
          - molr2*sth.hHsCoeffs_[coefLabel];

        iHsCoeffs_[coefLabel] =
            molr1*iHsCoeffs_[coefLabel]
          - molr2*sth.iHsCoeffs_[coefLabel];
     
        jHsCoeffs_[coefLabel] =
            molr1*jHsCoeffs_[coefLabel]
          - molr2*sth.jHsCoeffs_[coefLabel];

        kHsCoeffs_[coefLabel] =
            molr1*kHsCoeffs_[coefLabel]
          - molr2*sth.kHsCoeffs_[coefLabel];

        lHsCoeffs_[coefLabel] =
            molr1*lHsCoeffs_[coefLabel]
          - molr2*sth.lHsCoeffs_[coefLabel];
   
        mHsCoeffs_[coefLabel] =
            molr1*mHsCoeffs_[coefLabel]
          - molr2*sth.mHsCoeffs_[coefLabel];

        nHsCoeffs_[coefLabel] =
            molr1*nHsCoeffs_[coefLabel]
          - molr2*sth.nHsCoeffs_[coefLabel];
     
        oHsCoeffs_[coefLabel] =
            molr1*oHsCoeffs_[coefLabel]
          - molr2*sth.oHsCoeffs_[coefLabel];

        pHsCoeffs_[coefLabel] =
            molr1*pHsCoeffs_[coefLabel]
          - molr2*sth.pHsCoeffs_[coefLabel];
     
        qHsCoeffs_[coefLabel] =
            molr1*qHsCoeffs_[coefLabel]
          - molr2*sth.qHsCoeffs_[coefLabel];

        rHsCoeffs_[coefLabel] =
            molr1*rHsCoeffs_[coefLabel]
          - molr2*sth.rHsCoeffs_[coefLabel];


    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::splineThermo<EquationOfState> Foam::operator+
(
    const splineThermo<EquationOfState>& sth1,
    const splineThermo<EquationOfState>& sth2
)
{
    EquationOfState eofs = sth1;
    eofs += sth2;

    scalar molr1 = sth1.nMoles()/eofs.nMoles();
    scalar molr2 = sth2.nMoles()/eofs.nMoles();

    typename splineThermo<EquationOfState>::coeffArray a0CpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray aCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray bCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray cCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray dCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray eCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray fCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray gCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray hCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray iCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray jCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray kCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray lCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray mCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray nCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray oCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray pCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray qCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray rCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray a0HsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray aHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray bHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray cHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray dHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray eHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray fHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray gHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray hHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray iHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray jHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray kHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray lHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray mHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray nHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray oHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray pHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray qHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray rHsCoeffs;


    for
    (
        register label coefLabel=0;
        coefLabel<splineThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {

        a0CpCoeffs[coefLabel] =
            molr1*sth1.a0CpCoeffs_[coefLabel]
          + molr2*sth2.a0CpCoeffs_[coefLabel];

        aCpCoeffs[coefLabel] =
            molr1*sth1.aCpCoeffs_[coefLabel]
          + molr2*sth2.aCpCoeffs_[coefLabel];
        
        bCpCoeffs[coefLabel] =
            molr1*sth1.bCpCoeffs_[coefLabel]
          + molr2*sth2.bCpCoeffs_[coefLabel];

        cCpCoeffs[coefLabel] =
            molr1*sth1.cCpCoeffs_[coefLabel]
          + molr2*sth2.cCpCoeffs_[coefLabel];

        dCpCoeffs[coefLabel] =
            molr1*sth1.dCpCoeffs_[coefLabel]
          + molr2*sth2.dCpCoeffs_[coefLabel];

        eCpCoeffs[coefLabel] =
            molr1*sth1.eCpCoeffs_[coefLabel]
          + molr2*sth2.eCpCoeffs_[coefLabel];

        fCpCoeffs[coefLabel] =
            molr1*sth1.fCpCoeffs_[coefLabel]
          + molr2*sth2.fCpCoeffs_[coefLabel];
        
        gCpCoeffs[coefLabel] =
            molr1*sth1.gCpCoeffs_[coefLabel]
          + molr2*sth2.gCpCoeffs_[coefLabel];

        hCpCoeffs[coefLabel] =
            molr1*sth1.hCpCoeffs_[coefLabel]
          + molr2*sth2.hCpCoeffs_[coefLabel];

        iCpCoeffs[coefLabel] =
            molr1*sth1.iCpCoeffs_[coefLabel]
          + molr2*sth2.iCpCoeffs_[coefLabel];

        jCpCoeffs[coefLabel] =
            molr1*sth1.jCpCoeffs_[coefLabel]
          + molr2*sth2.jCpCoeffs_[coefLabel];

        kCpCoeffs[coefLabel] =
            molr1*sth1.kCpCoeffs_[coefLabel]
          + molr2*sth2.kCpCoeffs_[coefLabel];
        
        lCpCoeffs[coefLabel] =
            molr1*sth1.lCpCoeffs_[coefLabel]
          + molr2*sth2.lCpCoeffs_[coefLabel];

        mCpCoeffs[coefLabel] =
            molr1*sth1.mCpCoeffs_[coefLabel]
          + molr2*sth2.mCpCoeffs_[coefLabel];

        nCpCoeffs[coefLabel] =
            molr1*sth1.nCpCoeffs_[coefLabel]
          + molr2*sth2.nCpCoeffs_[coefLabel];

        oCpCoeffs[coefLabel] =
            molr1*sth1.oCpCoeffs_[coefLabel]
          + molr2*sth2.oCpCoeffs_[coefLabel];

        pCpCoeffs[coefLabel] =
            molr1*sth1.pCpCoeffs_[coefLabel]
          + molr2*sth2.pCpCoeffs_[coefLabel];

        qCpCoeffs[coefLabel] =
            molr1*sth1.qCpCoeffs_[coefLabel]
          + molr2*sth2.qCpCoeffs_[coefLabel];

        rCpCoeffs[coefLabel] =
            molr1*sth1.rCpCoeffs_[coefLabel]
          + molr2*sth2.rCpCoeffs_[coefLabel];

        a0HsCoeffs[coefLabel] =
            molr1*sth1.a0HsCoeffs_[coefLabel]
          + molr2*sth2.a0HsCoeffs_[coefLabel];

        aHsCoeffs[coefLabel] =
            molr1*sth1.aHsCoeffs_[coefLabel]
          + molr2*sth2.aHsCoeffs_[coefLabel];
        
        bHsCoeffs[coefLabel] =
            molr1*sth1.bHsCoeffs_[coefLabel]
          + molr2*sth2.bHsCoeffs_[coefLabel];

        cHsCoeffs[coefLabel] =
            molr1*sth1.cHsCoeffs_[coefLabel]
          + molr2*sth2.cHsCoeffs_[coefLabel];

        dHsCoeffs[coefLabel] =
            molr1*sth1.dHsCoeffs_[coefLabel]
          + molr2*sth2.dHsCoeffs_[coefLabel];

        eHsCoeffs[coefLabel] =
            molr1*sth1.eHsCoeffs_[coefLabel]
          + molr2*sth2.eHsCoeffs_[coefLabel];

        fHsCoeffs[coefLabel] =
            molr1*sth1.fHsCoeffs_[coefLabel]
          + molr2*sth2.fHsCoeffs_[coefLabel];
        
        gHsCoeffs[coefLabel] =
            molr1*sth1.gHsCoeffs_[coefLabel]
          + molr2*sth2.gHsCoeffs_[coefLabel];

        hHsCoeffs[coefLabel] =
            molr1*sth1.hHsCoeffs_[coefLabel]
          + molr2*sth2.hHsCoeffs_[coefLabel];

        iHsCoeffs[coefLabel] =
            molr1*sth1.iHsCoeffs_[coefLabel]
          + molr2*sth2.iHsCoeffs_[coefLabel];

        jHsCoeffs[coefLabel] =
            molr1*sth1.jHsCoeffs_[coefLabel]
          + molr2*sth2.jHsCoeffs_[coefLabel];

        kHsCoeffs[coefLabel] =
            molr1*sth1.kHsCoeffs_[coefLabel]
          + molr2*sth2.kHsCoeffs_[coefLabel];
        
        lHsCoeffs[coefLabel] =
            molr1*sth1.lHsCoeffs_[coefLabel]
          + molr2*sth2.lHsCoeffs_[coefLabel];

        mHsCoeffs[coefLabel] =
            molr1*sth1.mHsCoeffs_[coefLabel]
          + molr2*sth2.mHsCoeffs_[coefLabel];

        nHsCoeffs[coefLabel] =
            molr1*sth1.nHsCoeffs_[coefLabel]
          + molr2*sth2.nHsCoeffs_[coefLabel];

        oHsCoeffs[coefLabel] =
            molr1*sth1.oHsCoeffs_[coefLabel]
          + molr2*sth2.oHsCoeffs_[coefLabel];

        pHsCoeffs[coefLabel] =
            molr1*sth1.pHsCoeffs_[coefLabel]
          + molr2*sth2.pHsCoeffs_[coefLabel];

        qHsCoeffs[coefLabel] =
            molr1*sth1.qHsCoeffs_[coefLabel]
          + molr2*sth2.qHsCoeffs_[coefLabel];

        rHsCoeffs[coefLabel] =
            molr1*sth1.rHsCoeffs_[coefLabel]
          + molr2*sth2.rHsCoeffs_[coefLabel];
    }

    if
    (
        splineThermo<EquationOfState>::debug
     && notEqual(sth1.T0_, sth2.T0_)
    )
    {
        FatalErrorIn
        (
            "operator+"
            "(const splineThermo<EquationOfState>& sth1,"
            " const splineThermo<EquationOfState>& sth2)"
        )   << "T0 " << sth1.T0_ << " for "
            << (sth1.name().size() ? sth1.name() : "others")
            << " != " << sth2.T0_ << " for "
            << (sth2.name().size() ? sth2.name() : "others")
            << exit(FatalError);
    }

    return splineThermo<EquationOfState>
    (
        eofs,
        max(sth1.Tlow_, sth2.Tlow_),
        min(sth1.Thigh_, sth2.Thigh_),
        sth1.T0_,
        sth1.T1_,
        sth1.T2_,
        sth1.T3_,
        sth1.T4_,
        sth1.T5_,
        sth1.T6_,
        sth1.T7_,
        sth1.T8_,
        sth1.T9_,
        sth1.T10_,
        sth1.T11_,
        sth1.T12_,
        sth1.T13_,
        sth1.T14_,
        sth1.T15_,
        sth1.T16_,
        sth1.T17_,
        a0CpCoeffs,
        aCpCoeffs,
        bCpCoeffs,
        cCpCoeffs,
        dCpCoeffs,
        eCpCoeffs,
        fCpCoeffs,
        gCpCoeffs,
        hCpCoeffs,
        iCpCoeffs,
        jCpCoeffs,
        kCpCoeffs,
        lCpCoeffs,
        mCpCoeffs,
        nCpCoeffs,
        oCpCoeffs,
        pCpCoeffs,
        qCpCoeffs,
        rCpCoeffs,
        a0HsCoeffs,
        aHsCoeffs,
        bHsCoeffs,
        cHsCoeffs,
        dHsCoeffs,
        eHsCoeffs,
        fHsCoeffs,
        gHsCoeffs,
        hHsCoeffs,
        iHsCoeffs,
        jHsCoeffs,
        kHsCoeffs,
        lHsCoeffs,
        mHsCoeffs,
        nHsCoeffs,
        oHsCoeffs,
        pHsCoeffs,
        qHsCoeffs,
        rHsCoeffs
    );
}


template<class EquationOfState>
inline Foam::splineThermo<EquationOfState> Foam::operator-
(
    const splineThermo<EquationOfState>& sth1,
    const splineThermo<EquationOfState>& sth2
)
{
    EquationOfState eofs = sth1;
    eofs -= sth2;

    scalar molr1 = sth1.nMoles()/eofs.nMoles();
    scalar molr2 = sth2.nMoles()/eofs.nMoles();

    typename splineThermo<EquationOfState>::coeffArray a0CpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray aCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray bCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray cCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray dCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray eCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray fCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray gCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray hCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray iCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray jCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray kCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray lCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray mCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray nCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray oCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray pCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray qCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray rCpCoeffs;
    typename splineThermo<EquationOfState>::coeffArray a0HsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray aHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray bHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray cHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray dHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray eHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray fHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray gHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray hHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray iHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray jHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray kHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray lHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray mHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray nHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray oHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray pHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray qHsCoeffs;
    typename splineThermo<EquationOfState>::coeffArray rHsCoeffs;

    for
    (
        register label coefLabel=0;
        coefLabel<splineThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {

        a0CpCoeffs[coefLabel] =
            molr1*sth1.a0CpCoeffs_[coefLabel]
          - molr2*sth2.a0CpCoeffs_[coefLabel];

        aCpCoeffs[coefLabel] =
            molr1*sth1.aCpCoeffs_[coefLabel]
          - molr2*sth2.aCpCoeffs_[coefLabel];
 
        bCpCoeffs[coefLabel] =
            molr1*sth1.bCpCoeffs_[coefLabel]
          - molr2*sth2.bCpCoeffs_[coefLabel];

        cCpCoeffs[coefLabel] =
            molr1*sth1.cCpCoeffs_[coefLabel]
          - molr2*sth2.cCpCoeffs_[coefLabel];
              
        dCpCoeffs[coefLabel] =
            molr1*sth1.dCpCoeffs_[coefLabel]
          - molr2*sth2.dCpCoeffs_[coefLabel];

        eCpCoeffs[coefLabel] =
            molr1*sth1.eCpCoeffs_[coefLabel]
          - molr2*sth2.eCpCoeffs_[coefLabel];

        fCpCoeffs[coefLabel] =
            molr1*sth1.fCpCoeffs_[coefLabel]
          - molr2*sth2.fCpCoeffs_[coefLabel];
        
        gCpCoeffs[coefLabel] =
            molr1*sth1.gCpCoeffs_[coefLabel]
          - molr2*sth2.gCpCoeffs_[coefLabel];

        hCpCoeffs[coefLabel] =
            molr1*sth1.hCpCoeffs_[coefLabel]
          - molr2*sth2.hCpCoeffs_[coefLabel];

        iCpCoeffs[coefLabel] =
            molr1*sth1.iCpCoeffs_[coefLabel]
          - molr2*sth2.iCpCoeffs_[coefLabel];

        jCpCoeffs[coefLabel] =
            molr1*sth1.jCpCoeffs_[coefLabel]
          - molr2*sth2.jCpCoeffs_[coefLabel];

        kCpCoeffs[coefLabel] =
            molr1*sth1.kCpCoeffs_[coefLabel]
          - molr2*sth2.kCpCoeffs_[coefLabel];
        
        lCpCoeffs[coefLabel] =
            molr1*sth1.lCpCoeffs_[coefLabel]
          - molr2*sth2.lCpCoeffs_[coefLabel];

        mCpCoeffs[coefLabel] =
            molr1*sth1.mCpCoeffs_[coefLabel]
          - molr2*sth2.mCpCoeffs_[coefLabel];

        nCpCoeffs[coefLabel] =
            molr1*sth1.nCpCoeffs_[coefLabel]
          - molr2*sth2.nCpCoeffs_[coefLabel];

        oCpCoeffs[coefLabel] =
            molr1*sth1.oCpCoeffs_[coefLabel]
          - molr2*sth2.oCpCoeffs_[coefLabel];

        pCpCoeffs[coefLabel] =
            molr1*sth1.pCpCoeffs_[coefLabel]
          - molr2*sth2.pCpCoeffs_[coefLabel];

        qCpCoeffs[coefLabel] =
            molr1*sth1.qCpCoeffs_[coefLabel]
          - molr2*sth2.qCpCoeffs_[coefLabel];

        rCpCoeffs[coefLabel] =
            molr1*sth1.rCpCoeffs_[coefLabel]
          - molr2*sth2.rCpCoeffs_[coefLabel];

        a0HsCoeffs[coefLabel] =
            molr1*sth1.a0HsCoeffs_[coefLabel]
          - molr2*sth2.a0HsCoeffs_[coefLabel];

        aHsCoeffs[coefLabel] =
            molr1*sth1.aHsCoeffs_[coefLabel]
          - molr2*sth2.aHsCoeffs_[coefLabel];
        
        bHsCoeffs[coefLabel] =
            molr1*sth1.bHsCoeffs_[coefLabel]
          - molr2*sth2.bHsCoeffs_[coefLabel];

        cHsCoeffs[coefLabel] =
            molr1*sth1.cHsCoeffs_[coefLabel]
          - molr2*sth2.cHsCoeffs_[coefLabel];

        dHsCoeffs[coefLabel] =
            molr1*sth1.dHsCoeffs_[coefLabel]
          - molr2*sth2.dHsCoeffs_[coefLabel];

        eHsCoeffs[coefLabel] =
            molr1*sth1.eHsCoeffs_[coefLabel]
          - molr2*sth2.eHsCoeffs_[coefLabel];

        fHsCoeffs[coefLabel] =
            molr1*sth1.fHsCoeffs_[coefLabel]
          - molr2*sth2.fHsCoeffs_[coefLabel];
        
        gHsCoeffs[coefLabel] =
            molr1*sth1.gHsCoeffs_[coefLabel]
          - molr2*sth2.gHsCoeffs_[coefLabel];

        hHsCoeffs[coefLabel] =
            molr1*sth1.hHsCoeffs_[coefLabel]
          - molr2*sth2.hHsCoeffs_[coefLabel];

        iHsCoeffs[coefLabel] =
            molr1*sth1.iHsCoeffs_[coefLabel]
          - molr2*sth2.iHsCoeffs_[coefLabel];

        jHsCoeffs[coefLabel] =
            molr1*sth1.jHsCoeffs_[coefLabel]
          - molr2*sth2.jHsCoeffs_[coefLabel];

        kHsCoeffs[coefLabel] =
            molr1*sth1.kHsCoeffs_[coefLabel]
          - molr2*sth2.kHsCoeffs_[coefLabel];
        
        lHsCoeffs[coefLabel] =
            molr1*sth1.lHsCoeffs_[coefLabel]
          - molr2*sth2.lHsCoeffs_[coefLabel];

        mHsCoeffs[coefLabel] =
            molr1*sth1.mHsCoeffs_[coefLabel]
          - molr2*sth2.mHsCoeffs_[coefLabel];

        nHsCoeffs[coefLabel] =
            molr1*sth1.nHsCoeffs_[coefLabel]
          - molr2*sth2.nHsCoeffs_[coefLabel];

        oHsCoeffs[coefLabel] =
            molr1*sth1.oHsCoeffs_[coefLabel]
          - molr2*sth2.oHsCoeffs_[coefLabel];

        pHsCoeffs[coefLabel] =
            molr1*sth1.pHsCoeffs_[coefLabel]
          - molr2*sth2.pHsCoeffs_[coefLabel];

        qHsCoeffs[coefLabel] =
            molr1*sth1.qHsCoeffs_[coefLabel]
          - molr2*sth2.qHsCoeffs_[coefLabel];

        rHsCoeffs[coefLabel] =
            molr1*sth1.rHsCoeffs_[coefLabel]
          - molr2*sth2.rHsCoeffs_[coefLabel];
    }

    if
    (
        splineThermo<EquationOfState>::debug
     && notEqual(sth1.T0_, sth2.T0_)
    )
    {
        FatalErrorIn
        (
            "operator-"
            "(const splineThermo<EquationOfState>& st1,"
            " const splineThermo<EquationOfState>& st2)"
        )   << "T0 " << sth1.T0_ << " for "
            << (sth1.name().size() ? sth1.name() : "others")
            << " != " << sth2.T0_ << " for "
            << (sth2.name().size() ? sth2.name() : "others")
            << exit(FatalError);
    }

    return splineThermo<EquationOfState>
    (
        eofs,
        max(sth1.Tlow_, sth2.Tlow_),
        min(sth1.Thigh_, sth2.Thigh_),
        sth1.T0_,
        sth1.T1_,
        sth1.T2_,
        sth1.T3_,
        sth1.T4_,
        sth1.T5_,
        sth1.T6_,
        sth1.T7_,
        sth1.T8_,
        sth1.T9_,
        sth1.T10_,
        sth1.T11_,
        sth1.T12_,
        sth1.T13_,
        sth1.T14_,
        sth1.T15_,
        sth1.T16_,
        sth1.T17_,
        a0CpCoeffs,
        aCpCoeffs,
        bCpCoeffs,
        cCpCoeffs,
        dCpCoeffs,
        eCpCoeffs,
        fCpCoeffs,
        gCpCoeffs,
        hCpCoeffs,
        iCpCoeffs,
        jCpCoeffs,
        kCpCoeffs,
        lCpCoeffs,
        mCpCoeffs,
        nCpCoeffs,
        oCpCoeffs,
        pCpCoeffs,
        qCpCoeffs,
        rCpCoeffs,
        a0HsCoeffs,
        aHsCoeffs,
        bHsCoeffs,
        cHsCoeffs,
        dHsCoeffs,
        eHsCoeffs,
        fHsCoeffs,
        gHsCoeffs,
        hHsCoeffs,
        iHsCoeffs,
        jHsCoeffs,
        kHsCoeffs,
        lHsCoeffs,
        mHsCoeffs,
        nHsCoeffs,
        oHsCoeffs,
        pHsCoeffs,
        qHsCoeffs,
        rHsCoeffs
    );
}


template<class EquationOfState>
inline Foam::splineThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const splineThermo<EquationOfState>& sth
)
{
    return splineThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(sth),
        sth.Tlow_,
        sth.Thigh_,
        sth.T0_,
        sth.T1_,
        sth.T2_,
        sth.T3_,
        sth.T4_,
        sth.T5_,
        sth.T6_,
        sth.T7_,
        sth.T8_,
        sth.T9_,
        sth.T10_,
        sth.T11_,
        sth.T12_,
        sth.T13_,
        sth.T14_,
        sth.T15_,
        sth.T16_,
        sth.T17_,
        sth.a0CpCoeffs_,
        sth.aCpCoeffs_,
        sth.bCpCoeffs_,
        sth.cCpCoeffs_,
        sth.dCpCoeffs_,
        sth.eCpCoeffs_,
        sth.fCpCoeffs_,
        sth.gCpCoeffs_,
        sth.hCpCoeffs_,
        sth.iCpCoeffs_,
        sth.jCpCoeffs_,
        sth.kCpCoeffs_,
        sth.lCpCoeffs_,
        sth.mCpCoeffs_,
        sth.nCpCoeffs_,
        sth.oCpCoeffs_,
        sth.pCpCoeffs_,
        sth.qCpCoeffs_,
        sth.rCpCoeffs_,
        sth.a0HsCoeffs_,
        sth.aHsCoeffs_,
        sth.bHsCoeffs_,
        sth.cHsCoeffs_,
        sth.dHsCoeffs_,
        sth.eHsCoeffs_,
        sth.fHsCoeffs_,
        sth.gHsCoeffs_,
        sth.hHsCoeffs_,
        sth.iHsCoeffs_,
        sth.jHsCoeffs_,
        sth.kHsCoeffs_,
        sth.lHsCoeffs_,
        sth.mHsCoeffs_,
        sth.nHsCoeffs_,
        sth.oHsCoeffs_,
        sth.pHsCoeffs_,
        sth.qHsCoeffs_,
        sth.rHsCoeffs_
    );
}


template<class EquationOfState>
inline Foam::splineThermo<EquationOfState> Foam::operator==
(
    const splineThermo<EquationOfState>& sth1,
    const splineThermo<EquationOfState>& sth2
)
{
    return sth2 - sth1;
}


// ************************************************************************* //
