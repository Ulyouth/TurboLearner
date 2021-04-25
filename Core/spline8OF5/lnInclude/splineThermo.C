/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
void Foam::splineThermo<EquationOfState>::checkInputData() const
{
    if (Tlow_ >= Thigh_)
    {
        FatalErrorIn("splineThermo<EquationOfState>::check()")
            << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }

    if (T1_ <= Tlow_)
    {
        FatalErrorIn("splineThermo<EquationOfState>::check()")
            << "T1(" << T1_ << ") <= Tlow(" << Tlow_ << ')'
            << exit(FatalError);
    }

    if (T16_ > Thigh_)
    {
        FatalErrorIn("splineThermo<EquationOfState>::check()")
            << "T16(" << T16_ << ") > Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::splineThermo<EquationOfState>::splineThermo(Istream& is)
:
    EquationOfState(is),
    Tlow_(readScalar(is)),
    Thigh_(readScalar(is)),
    T0_(readScalar(is)),
    T1_(readScalar(is)),
    T2_(readScalar(is)),
    T3_(readScalar(is)),
    T4_(readScalar(is)),
    T5_(readScalar(is)),
    T6_(readScalar(is)),
    T7_(readScalar(is)),
    T8_(readScalar(is)),
    T9_(readScalar(is)),
    T10_(readScalar(is)),
    T11_(readScalar(is)),
    T12_(readScalar(is)),
    T13_(readScalar(is)),
    T14_(readScalar(is)),
    T15_(readScalar(is)),
    T16_(readScalar(is)),
    T17_(readScalar(is))

{
    checkInputData();

    forAll(a0CpCoeffs_, i)
    {
        is >> a0CpCoeffs_[i];
    }

    forAll(aCpCoeffs_, i)
    {
        is >> aCpCoeffs_[i];
    }

    forAll(bCpCoeffs_, i)
    {
        is >> bCpCoeffs_[i];
    }

    forAll(cCpCoeffs_, i)
    {
        is >> cCpCoeffs_[i];
    }

    forAll(dCpCoeffs_, i)
    {
        is >> dCpCoeffs_[i];
    }
    
    forAll(eCpCoeffs_, i)
    {
        is >> eCpCoeffs_[i];
    }

    forAll(fCpCoeffs_, i)
    {
        is >> fCpCoeffs_[i];
    }

    forAll(gCpCoeffs_, i)
    {
        is >> gCpCoeffs_[i];
    }

    forAll(hCpCoeffs_, i)
    {
        is >> hCpCoeffs_[i];
    }

    forAll(iCpCoeffs_, i)
    {
        is >> iCpCoeffs_[i];
    }

    forAll(jCpCoeffs_, i)
    {
        is >> jCpCoeffs_[i];
    }
    
    forAll(kCpCoeffs_, i)
    {
        is >> kCpCoeffs_[i];
    }

    forAll(lCpCoeffs_, i)
    {
        is >> lCpCoeffs_[i];
    }

    forAll(mCpCoeffs_, i)
    {
        is >> mCpCoeffs_[i];
    }

    forAll(nCpCoeffs_, i)
    {
        is >> nCpCoeffs_[i];
    }

    forAll(oCpCoeffs_, i)
    {
        is >> oCpCoeffs_[i];
    }
    
    forAll(pCpCoeffs_, i)
    {
        is >> pCpCoeffs_[i];
    }

    forAll(qCpCoeffs_, i)
    {
        is >> qCpCoeffs_[i];
    }

    forAll(rCpCoeffs_, i)
    {
        is >> rCpCoeffs_[i];
    }

   forAll(a0HsCoeffs_, i)
    {
        is >> a0HsCoeffs_[i];
    }

   forAll(aHsCoeffs_, i)
    {
        is >> aHsCoeffs_[i];
    }

    forAll(bHsCoeffs_, i)
    {
        is >> bHsCoeffs_[i];
    }

    forAll(cHsCoeffs_, i)
    {
        is >> cHsCoeffs_[i];
    }

    forAll(dHsCoeffs_, i)
    {
        is >> dHsCoeffs_[i];
    }
    
    forAll(eHsCoeffs_, i)
    {
        is >> eHsCoeffs_[i];
    }

    forAll(fHsCoeffs_, i)
    {
        is >> fHsCoeffs_[i];
    }

    forAll(gHsCoeffs_, i)
    {
        is >> gHsCoeffs_[i];
    }

    forAll(hHsCoeffs_, i)
    {
        is >> hHsCoeffs_[i];
    }

    forAll(iHsCoeffs_, i)
    {
        is >> iHsCoeffs_[i];
    }

    forAll(jHsCoeffs_, i)
    {
        is >> jHsCoeffs_[i];
    }
    
    forAll(kHsCoeffs_, i)
    {
        is >> kHsCoeffs_[i];
    }

    forAll(lHsCoeffs_, i)
    {
        is >> lHsCoeffs_[i];
    }

    forAll(mHsCoeffs_, i)
    {
        is >> mHsCoeffs_[i];
    }

    forAll(nHsCoeffs_, i)
    {
        is >> nHsCoeffs_[i];
    }

    forAll(oHsCoeffs_, i)
    {
        is >> oHsCoeffs_[i];
    }
    
    forAll(pHsCoeffs_, i)
    {
        is >> pHsCoeffs_[i];
    }

    forAll(qHsCoeffs_, i)
    {
        is >> qHsCoeffs_[i];
    }

    forAll(rHsCoeffs_, i)
    {
        is >> rHsCoeffs_[i];
    }
    // Check state of Istream
    is.check("splineThermo::splineThermo(Istream& is)");
}


template<class EquationOfState>
Foam::splineThermo<EquationOfState>::splineThermo(const dictionary& dict)
:
    EquationOfState(dict),
    Tlow_(readScalar(dict.subDict("thermodynamics").lookup("Tlow"))),
    Thigh_(readScalar(dict.subDict("thermodynamics").lookup("Thigh"))),
    T0_(readScalar(dict.subDict("thermodynamics").lookup("T0"))),
    T1_(readScalar(dict.subDict("thermodynamics").lookup("T1"))),
    T2_(readScalar(dict.subDict("thermodynamics").lookup("T2"))),
    T3_(readScalar(dict.subDict("thermodynamics").lookup("T3"))),
    T4_(readScalar(dict.subDict("thermodynamics").lookup("T4"))),
    T5_(readScalar(dict.subDict("thermodynamics").lookup("T5"))),
    T6_(readScalar(dict.subDict("thermodynamics").lookup("T6"))),
    T7_(readScalar(dict.subDict("thermodynamics").lookup("T7"))),
    T8_(readScalar(dict.subDict("thermodynamics").lookup("T8"))),
    T9_(readScalar(dict.subDict("thermodynamics").lookup("T9"))),
    T10_(readScalar(dict.subDict("thermodynamics").lookup("T10"))),
    T11_(readScalar(dict.subDict("thermodynamics").lookup("T11"))),
    T12_(readScalar(dict.subDict("thermodynamics").lookup("T12"))),
    T13_(readScalar(dict.subDict("thermodynamics").lookup("T13"))),
    T14_(readScalar(dict.subDict("thermodynamics").lookup("T14"))),
    T15_(readScalar(dict.subDict("thermodynamics").lookup("T15"))),
    T16_(readScalar(dict.subDict("thermodynamics").lookup("T16"))),
    a0CpCoeffs_(dict.subDict("thermodynamics").lookup("a0CpCoeffs")),
    aCpCoeffs_(dict.subDict("thermodynamics").lookup("aCpCoeffs")),
    bCpCoeffs_(dict.subDict("thermodynamics").lookup("bCpCoeffs")),
    cCpCoeffs_(dict.subDict("thermodynamics").lookup("cCpCoeffs")),
    dCpCoeffs_(dict.subDict("thermodynamics").lookup("dCpCoeffs")),
    eCpCoeffs_(dict.subDict("thermodynamics").lookup("eCpCoeffs")),
    fCpCoeffs_(dict.subDict("thermodynamics").lookup("fCpCoeffs")),
    gCpCoeffs_(dict.subDict("thermodynamics").lookup("gCpCoeffs")),
    hCpCoeffs_(dict.subDict("thermodynamics").lookup("hCpCoeffs")),
    iCpCoeffs_(dict.subDict("thermodynamics").lookup("iCpCoeffs")),
    jCpCoeffs_(dict.subDict("thermodynamics").lookup("jCpCoeffs")),
    kCpCoeffs_(dict.subDict("thermodynamics").lookup("kCpCoeffs")),
    lCpCoeffs_(dict.subDict("thermodynamics").lookup("lCpCoeffs")),
    mCpCoeffs_(dict.subDict("thermodynamics").lookup("mCpCoeffs")),
    nCpCoeffs_(dict.subDict("thermodynamics").lookup("nCpCoeffs")),
    oCpCoeffs_(dict.subDict("thermodynamics").lookup("oCpCoeffs")),
    pCpCoeffs_(dict.subDict("thermodynamics").lookup("pCpCoeffs")),
    qCpCoeffs_(dict.subDict("thermodynamics").lookup("qCpCoeffs")),
    rCpCoeffs_(dict.subDict("thermodynamics").lookup("rCpCoeffs")),
    a0HsCoeffs_(dict.subDict("thermodynamics").lookup("a0hsCoeffs")),
    aHsCoeffs_(dict.subDict("thermodynamics").lookup("ahsCoeffs")),
    bHsCoeffs_(dict.subDict("thermodynamics").lookup("bhsCoeffs")),
    cHsCoeffs_(dict.subDict("thermodynamics").lookup("chsCoeffs")),
    dHsCoeffs_(dict.subDict("thermodynamics").lookup("dhsCoeffs")),
    eHsCoeffs_(dict.subDict("thermodynamics").lookup("ehsCoeffs")),
    fHsCoeffs_(dict.subDict("thermodynamics").lookup("fhsCoeffs")),
    gHsCoeffs_(dict.subDict("thermodynamics").lookup("ghsCoeffs")),
    hHsCoeffs_(dict.subDict("thermodynamics").lookup("hhsCoeffs")),
    iHsCoeffs_(dict.subDict("thermodynamics").lookup("ihsCoeffs")),
    jHsCoeffs_(dict.subDict("thermodynamics").lookup("jhsCoeffs")),
    kHsCoeffs_(dict.subDict("thermodynamics").lookup("khsCoeffs")),
    lHsCoeffs_(dict.subDict("thermodynamics").lookup("lhsCoeffs")),
    mHsCoeffs_(dict.subDict("thermodynamics").lookup("mhsCoeffs")),
    nHsCoeffs_(dict.subDict("thermodynamics").lookup("nhsCoeffs")),
    oHsCoeffs_(dict.subDict("thermodynamics").lookup("ohsCoeffs")),
    pHsCoeffs_(dict.subDict("thermodynamics").lookup("phsCoeffs")),
    qHsCoeffs_(dict.subDict("thermodynamics").lookup("qhsCoeffs")),
    rHsCoeffs_(dict.subDict("thermodynamics").lookup("rhsCoeffs"))
{
    checkInputData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::splineThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("Tlow", Tlow_);
    dict.add("Thigh", Thigh_);
    dict.add("T0", T0_);
    dict.add("T1", T1_);
    dict.add("T2", T2_);
    dict.add("T3", T3_);
    dict.add("T4", T4_);
    dict.add("T5", T5_);
    dict.add("T6", T6_);
    dict.add("T7", T7_);
    dict.add("T8", T8_);
    dict.add("T9", T9_);
    dict.add("T10", T10_);
    dict.add("T11", T11_);
    dict.add("T12", T12_);
    dict.add("T13", T13_);
    dict.add("T14", T14_);
    dict.add("T15", T15_);
    dict.add("T16", T16_);
    dict.add("T17", T17_);
    dict.add("a0CpCoeffs", a0CpCoeffs_);
    dict.add("aCpCoeffs", aCpCoeffs_);
    dict.add("bCpCoeffs", bCpCoeffs_);
    dict.add("cCpCoeffs", cCpCoeffs_);
    dict.add("dCpCoeffs", dCpCoeffs_);
    dict.add("eCpCoeffs", eCpCoeffs_);
    dict.add("fCpCoeffs", fCpCoeffs_);
    dict.add("gCpCoeffs", gCpCoeffs_);
    dict.add("hCpCoeffs", hCpCoeffs_);
    dict.add("iCpCoeffs", iCpCoeffs_);
    dict.add("jCpCoeffs", jCpCoeffs_);
    dict.add("kCpCoeffs", kCpCoeffs_);
    dict.add("lCpCoeffs", lCpCoeffs_);
    dict.add("mCpCoeffs", mCpCoeffs_);
    dict.add("nCpCoeffs", nCpCoeffs_);
    dict.add("oCpCoeffs", oCpCoeffs_);
    dict.add("pCpCoeffs", pCpCoeffs_);
    dict.add("qCpCoeffs", qCpCoeffs_);
    dict.add("rCpCoeffs", rCpCoeffs_);
    dict.add("a0hsCoeffs", a0HsCoeffs_);
    dict.add("ahsCoeffs", aHsCoeffs_);
    dict.add("bhsCoeffs", bHsCoeffs_);
    dict.add("chsCoeffs", cHsCoeffs_);
    dict.add("dhsCoeffs", dHsCoeffs_);
    dict.add("ehsCoeffs", eHsCoeffs_);
    dict.add("fhsCoeffs", fHsCoeffs_);
    dict.add("ghsCoeffs", gHsCoeffs_);
    dict.add("hhsCoeffs", hHsCoeffs_);
    dict.add("ihsCoeffs", iHsCoeffs_);
    dict.add("jhsCoeffs", jHsCoeffs_);
    dict.add("khsCoeffs", kHsCoeffs_);
    dict.add("lhsCoeffs", lHsCoeffs_);
    dict.add("mhsCoeffs", mHsCoeffs_);
    dict.add("nhsCoeffs", nHsCoeffs_);
    dict.add("ohsCoeffs", oHsCoeffs_);
    dict.add("phsCoeffs", pHsCoeffs_);
    dict.add("qhsCoeffs", qHsCoeffs_);
    dict.add("rhsCoeffs", rHsCoeffs_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const splineThermo<EquationOfState>& sth
)
{
    os  << static_cast<const EquationOfState&>(sth) << nl
        << "    " << sth.Tlow_
        << tab << sth.Thigh_
        << tab << sth.T0_
        << tab << sth.T1_
        << tab << sth.T2_
        << tab << sth.T3_
        << tab << sth.T4_
        << tab << sth.T5_
        << tab << sth.T6_
        << tab << sth.T7_
        << tab << sth.T8_
        << tab << sth.T9_
        << tab << sth.T10_
        << tab << sth.T11_
        << tab << sth.T12_
        << tab << sth.T13_
        << tab << sth.T14_
        << tab << sth.T15_
        << tab << sth.T16_
        << tab << sth.T17_;

    os << nl << "    ";

    forAll(sth.a0CpCoeffs_, i)
    {
        os << sth.a0CpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.aCpCoeffs_, i)
    {
        os << sth.aCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.bCpCoeffs_, i)
    {
        os << sth.bCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.cCpCoeffs_, i)
    {
        os << sth.cCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.dCpCoeffs_, i)
    {
        os << sth.dCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.eCpCoeffs_, i)
    {
        os << sth.eCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.fCpCoeffs_, i)
    {
        os << sth.fCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.gCpCoeffs_, i)
    {
        os << sth.gCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.hCpCoeffs_, i)
    {
        os << sth.hCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.iCpCoeffs_, i)
    {
        os << sth.iCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.jCpCoeffs_, i)
    {
        os << sth.jCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.kCpCoeffs_, i)
    {
        os << sth.kCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.lCpCoeffs_, i)
    {
        os << sth.lCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.mCpCoeffs_, i)
    {
        os << sth.mCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.nCpCoeffs_, i)
    {
        os << sth.nCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.oCpCoeffs_, i)
    {
        os << sth.oCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.pCpCoeffs_, i)
    {
        os << sth.pCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.qCpCoeffs_, i)
    {
        os << sth.qCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.rCpCoeffs_, i)
    {
        os << sth.rCpCoeffs_[i] << ' ';
    }


    os << nl << "    ";

    forAll(sth.a0hsCoeffs_, i)
    {
        os << sth.a0hsCoeffs_[i] << ' ';
    }


    os << nl << "    ";

    forAll(sth.ahsCoeffs_, i)
    {
        os << sth.ahsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.bhsCoeffs_, i)
    {
        os << sth.bhsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.chsCoeffs_, i)
    {
        os << sth.chsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.dhsCoeffs_, i)
    {
        os << sth.dhsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.ehsCoeffs_, i)
    {
        os << sth.ehsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.fhsCoeffs_, i)
    {
        os << sth.fhsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.ghsCoeffs_, i)
    {
        os << sth.ghsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.hhsCoeffs_, i)
    {
        os << sth.hhsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.ihsCoeffs_, i)
    {
        os << sth.ihsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.jhsCoeffs_, i)
    {
        os << sth.jhsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.khsCoeffs_, i)
    {
        os << sth.khsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.lhsCoeffs_, i)
    {
        os << sth.lhsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.mhsCoeffs_, i)
    {
        os << sth.mhsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.nhsCoeffs_, i)
    {
        os << sth.nhsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.ohsCoeffs_, i)
    {
        os << sth.ohsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.phsCoeffs_, i)
    {
        os << sth.phsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.qhsCoeffs_, i)
    {
        os << sth.qhsCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(sth.rhsCoeffs_, i)
    {
        os << sth.rhsCoeffs_[i] << ' ';
    }


    os << endl;

    os.check
    (
        "operator<<(Ostream& os, const splineThermo<EquationOfState>& sth)"
    );

    return os;
}


// ************************************************************************* //
