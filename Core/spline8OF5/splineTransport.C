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
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
void Foam::splineTransport<Thermo>::checkInputData() const
{
    if (Tlow_ >= Thigh_)
    {
        FatalErrorIn("splineTransport<Thermo>::check()")
            << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }

    if (T1_ <= Tlow_)
    {
        FatalErrorIn("splineTransport<Thermo>::check()")
            << "T1(" << T1_ << ") <= Tlow(" << Tlow_ << ')'
            << exit(FatalError);
    }

    if (T15_ > Thigh_)
    {
        FatalErrorIn("splineTransport<Thermo>::check()")
            << "T15(" << T15_ << ") > Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::splineTransport<Thermo>::splineTransport(Istream& is)
:
    Thermo(is),
    Tlow_(readScalar(is)),
    Thigh_(readScalar(is)),
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
    T15_(readScalar(is))
    
{
    checkInputData();

    forAll(amuCoeffs_, i)
    {
        is >> amuCoeffs_[i];
    }

    forAll(bmuCoeffs_, i)
    {
        is >> bmuCoeffs_[i];
    }

    forAll(cmuCoeffs_, i)
    {
        is >> cmuCoeffs_[i];
    }

    forAll(dmuCoeffs_, i)
    {
        is >> dmuCoeffs_[i];
    }
    forAll(emuCoeffs_, i)
    {
        is >> emuCoeffs_[i];
    }
    forAll(fmuCoeffs_, i)
    {
        is >> fmuCoeffs_[i];
    }

    forAll(gmuCoeffs_, i)
    {
        is >> gmuCoeffs_[i];
    }

    forAll(hmuCoeffs_, i)
    {
        is >> hmuCoeffs_[i];
    }

    forAll(imuCoeffs_, i)
    {
        is >> imuCoeffs_[i];
    }
    forAll(jmuCoeffs_, i)
    {
        is >> jmuCoeffs_[i];
    }
    forAll(kmuCoeffs_, i)
    {
        is >> kmuCoeffs_[i];
    }

    forAll(lmuCoeffs_, i)
    {
        is >> lmuCoeffs_[i];
    }

    forAll(mmuCoeffs_, i)
    {
        is >> mmuCoeffs_[i];
    }

    forAll(nmuCoeffs_, i)
    {
        is >> nmuCoeffs_[i];
    }

    forAll(omuCoeffs_, i)
    {
        is >> omuCoeffs_[i];
    }

    forAll(pmuCoeffs_, i)
    {
        is >> pmuCoeffs_[i];
    }

    forAll(akappaCoeffs_, i)
    {
        is >> akappaCoeffs_[i];
    }

    forAll(bkappaCoeffs_, i)
    {
        is >> bkappaCoeffs_[i];
    }

    forAll(ckappaCoeffs_, i)
    {
        is >> ckappaCoeffs_[i];
    }

    forAll(dkappaCoeffs_, i)
    {
        is >> dkappaCoeffs_[i];
    }
    forAll(ekappaCoeffs_, i)
    {
        is >> ekappaCoeffs_[i];
    }
    forAll(fkappaCoeffs_, i)
    {
        is >> fkappaCoeffs_[i];
    }

    forAll(gkappaCoeffs_, i)
    {
        is >> gkappaCoeffs_[i];
    }

    forAll(hkappaCoeffs_, i)
    {
        is >> hkappaCoeffs_[i];
    }

    forAll(ikappaCoeffs_, i)
    {
        is >> ikappaCoeffs_[i];
    }
    forAll(jkappaCoeffs_, i)
    {
        is >> jkappaCoeffs_[i];
    }
    forAll(kkappaCoeffs_, i)
    {
        is >> kkappaCoeffs_[i];
    }

    forAll(lkappaCoeffs_, i)
    {
        is >> lkappaCoeffs_[i];
    }

    forAll(mkappaCoeffs_, i)
    {
        is >> mkappaCoeffs_[i];
    }

    forAll(nkappaCoeffs_, i)
    {
        is >> nkappaCoeffs_[i];
    }
    forAll(okappaCoeffs_, i)
    {
        is >> okappaCoeffs_[i];
    }
    forAll(pkappaCoeffs_, i)
    {
        is >> pkappaCoeffs_[i];
    }
    // Check state of Istream
    is.check("splineTransport::splineTransport(Istream& is)");
}

template<class Thermo>
Foam::splineTransport<Thermo>::splineTransport(const dictionary& dict)
:
    Thermo(dict),
    Tlow_(readScalar(dict.subDict("transport").lookup("Tlow"))),
    Thigh_(readScalar(dict.subDict("transport").lookup("Thigh"))),
    T1_(readScalar(dict.subDict("transport").lookup("T1"))),
    T2_(readScalar(dict.subDict("transport").lookup("T2"))),
    T3_(readScalar(dict.subDict("transport").lookup("T3"))),
    T4_(readScalar(dict.subDict("transport").lookup("T4"))),
    T5_(readScalar(dict.subDict("transport").lookup("T5"))),
    T6_(readScalar(dict.subDict("transport").lookup("T6"))),
    T7_(readScalar(dict.subDict("transport").lookup("T7"))),
    T8_(readScalar(dict.subDict("transport").lookup("T8"))),
    T9_(readScalar(dict.subDict("transport").lookup("T9"))),
    T10_(readScalar(dict.subDict("transport").lookup("T10"))),
    T11_(readScalar(dict.subDict("transport").lookup("T11"))),
    T12_(readScalar(dict.subDict("transport").lookup("T12"))),
    T13_(readScalar(dict.subDict("transport").lookup("T13"))),
    T14_(readScalar(dict.subDict("transport").lookup("T14"))),
    T15_(readScalar(dict.subDict("transport").lookup("T15"))),
    amuCoeffs_(dict.subDict("transport").lookup("amuCoeffs")),
    bmuCoeffs_(dict.subDict("transport").lookup("bmuCoeffs")),
    cmuCoeffs_(dict.subDict("transport").lookup("cmuCoeffs")),
    dmuCoeffs_(dict.subDict("transport").lookup("dmuCoeffs")),
    emuCoeffs_(dict.subDict("transport").lookup("emuCoeffs")),
    fmuCoeffs_(dict.subDict("transport").lookup("fmuCoeffs")),
    gmuCoeffs_(dict.subDict("transport").lookup("gmuCoeffs")),
    hmuCoeffs_(dict.subDict("transport").lookup("hmuCoeffs")),
    imuCoeffs_(dict.subDict("transport").lookup("imuCoeffs")),
    jmuCoeffs_(dict.subDict("transport").lookup("jmuCoeffs")),
    kmuCoeffs_(dict.subDict("transport").lookup("kmuCoeffs")),
    lmuCoeffs_(dict.subDict("transport").lookup("lmuCoeffs")),
    mmuCoeffs_(dict.subDict("transport").lookup("mmuCoeffs")),
    nmuCoeffs_(dict.subDict("transport").lookup("nmuCoeffs")),
    omuCoeffs_(dict.subDict("transport").lookup("omuCoeffs")),
    pmuCoeffs_(dict.subDict("transport").lookup("pmuCoeffs")),
    akappaCoeffs_(dict.subDict("transport").lookup("akappaCoeffs")),
    bkappaCoeffs_(dict.subDict("transport").lookup("bkappaCoeffs")),
    ckappaCoeffs_(dict.subDict("transport").lookup("ckappaCoeffs")),
    dkappaCoeffs_(dict.subDict("transport").lookup("dkappaCoeffs")),
    ekappaCoeffs_(dict.subDict("transport").lookup("ekappaCoeffs")),
    fkappaCoeffs_(dict.subDict("transport").lookup("fkappaCoeffs")),
    gkappaCoeffs_(dict.subDict("transport").lookup("gkappaCoeffs")),
    hkappaCoeffs_(dict.subDict("transport").lookup("hkappaCoeffs")),
    ikappaCoeffs_(dict.subDict("transport").lookup("ikappaCoeffs")),
    jkappaCoeffs_(dict.subDict("transport").lookup("jkappaCoeffs")),
    kkappaCoeffs_(dict.subDict("transport").lookup("kkappaCoeffs")),
    lkappaCoeffs_(dict.subDict("transport").lookup("lkappaCoeffs")),
    mkappaCoeffs_(dict.subDict("transport").lookup("mkappaCoeffs")),
    nkappaCoeffs_(dict.subDict("transport").lookup("nkappaCoeffs")),
    okappaCoeffs_(dict.subDict("transport").lookup("okappaCoeffs")),
    pkappaCoeffs_(dict.subDict("transport").lookup("pkappaCoeffs"))
 

{
    checkInputData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::splineTransport<Thermo>::write(Ostream& os) const
{
    Thermo::write(os);

    dictionary dict("transport");
    dict.add("Tlow", Tlow_);
    dict.add("Thigh", Thigh_);
    dict.add("T1", T1_);
    dict.add("T2", T2_);
    dict.add("T3", T3_);
    dict.add("T4", T4_);
    dict.add("T5", T1_);
    dict.add("T6", T2_);
    dict.add("T7", T3_);
    dict.add("T8", T4_);
    dict.add("T9", T1_);
    dict.add("T10", T2_);
    dict.add("T11", T3_);
    dict.add("T12", T4_);
    dict.add("T13", T1_);
    dict.add("T14", T2_);
    dict.add("T15", T3_);
    dict.add("amuCoeffs", amuCoeffs_);
    dict.add("bmuCoeffs", bmuCoeffs_);
    dict.add("cmuCoeffs", cmuCoeffs_);
    dict.add("dmuCoeffs", dmuCoeffs_);
    dict.add("emuCoeffs", emuCoeffs_);
    dict.add("fmuCoeffs", fmuCoeffs_);
    dict.add("gmuCoeffs", gmuCoeffs_);
    dict.add("hmuCoeffs", hmuCoeffs_);
    dict.add("imuCoeffs", imuCoeffs_);
    dict.add("jmuCoeffs", jmuCoeffs_);
    dict.add("kmuCoeffs", kmuCoeffs_);
    dict.add("lmuCoeffs", lmuCoeffs_);
    dict.add("mmuCoeffs", mmuCoeffs_);
    dict.add("nmuCoeffs", nmuCoeffs_);
    dict.add("omuCoeffs", omuCoeffs_);
    dict.add("pmuCoeffs", pmuCoeffs_);
    dict.add("akappaCoeffs", akappaCoeffs_);
    dict.add("bkappaCoeffs", bkappaCoeffs_);
    dict.add("ckappaCoeffs", ckappaCoeffs_);
    dict.add("dkappaCoeffs", dkappaCoeffs_);
    dict.add("ekappaCoeffs", ekappaCoeffs_);
    dict.add("fkappaCoeffs", fkappaCoeffs_);
    dict.add("gkappaCoeffs", gkappaCoeffs_);
    dict.add("hkappaCoeffs", hkappaCoeffs_);
    dict.add("ikappaCoeffs", ikappaCoeffs_);
    dict.add("jkappaCoeffs", jkappaCoeffs_);
    dict.add("kkappaCoeffs", kkappaCoeffs_);
    dict.add("lkappaCoeffs", lkappaCoeffs_);
    dict.add("mkappaCoeffs", mkappaCoeffs_);
    dict.add("nkappaCoeffs", nkappaCoeffs_);
    dict.add("okappaCoeffs", okappaCoeffs_);
    dict.add("pkappaCoeffs", pkappaCoeffs_);


    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * IOstream Operator  * * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<

(Ostream& os, const splineTransport<Thermo>& spt)
        
{
    os  << static_cast<const Thermo&>(spt) << nl
        << "    " << spt.Tlow_
        << tab << spt.Thigh_
        << tab << spt.T1_
        << tab << spt.T2_
        << tab << spt.T3_
        << tab << spt.T4_
        << tab << spt.T5_
        << tab << spt.T6_
        << tab << spt.T7_
        << tab << spt.T8_
        << tab << spt.T9_
        << tab << spt.T10_
        << tab << spt.T11_
        << tab << spt.T12_
        << tab << spt.T13_
        << tab << spt.T14_
        << tab << spt.T15_;

    os << nl << "    ";

    forAll(spt.amuCoeffs_, i)
    {
        os << spt.amuCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.bmuCoeffs_, i)
    {
        os << spt.bmuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.cmuCoeffs_, i)
    {
        os << spt.cmuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.dmuCoeffs_, i)
    {
        os << spt.dmuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.emuCoeffs_, i)
    {
        os << spt.emuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.fmuCoeffs_, i)
    {
        os << spt.fmuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.gmuCoeffs_, i)
    {
        os << spt.gmuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.hmuCoeffs_, i)
    {
        os << spt.hmuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.imuCoeffs_, i)
    {
        os << spt.imuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.jmuCoeffs_, i)
    {
        os << spt.jmuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.kmuCoeffs_, i)
    {
        os << spt.kmuCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.lmuCoeffs_, i)
    {
        os << spt.lmuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.mmuCoeffs_, i)
    {
        os << spt.mmuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.nmuCoeffs_, i)
    {
        os << spt.nmuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.omuCoeffs_, i)
    {
        os << spt.omuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.pmuCoeffs_, i)
    {
        os << spt.pmuCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.akappaCoeffs_, i)
    {
        os << spt.akappaCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.bkappaCoeffs_, i)
    {
        os << spt.bkappaCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.ckappaCoeffs_, i)
    {
        os << spt.ckappaCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.dkappaCoeffs_, i)
    {
        os << spt.dkappaCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.ekappaCoeffs_, i)
    {
        os << spt.ekappaCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.fkappaCoeffs_, i)
    {
        os << spt.fkappaCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.gkappaCoeffs_, i)
    {
        os << spt.gkappaCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.hkappaCoeffs_, i)
    {
        os << spt.hkappaCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.ikappaCoeffs_, i)
    {
        os << spt.ikappaCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.jkappaCoeffs_, i)
    {
        os << spt.jkappaCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.kkappaCoeffs_, i)
    {
        os << spt.kkappaCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.lkappaCoeffs_, i)
    {
        os << spt.lkappaCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.mkappaCoeffs_, i)
    {
        os << spt.mkappaCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.nkappaCoeffs_, i)
    {
        os << spt.nkappaCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.okappaCoeffs_, i)
    {
        os << spt.okappaCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(spt.pkappaCoeffs_, i)
    {
        os << spt.pkappaCoeffs_[i] << ' ';
    }

    os.check
    (
        "Ostream& operator<<(Ostream& os, const splineTransport<Thermo>& spt)"
    );

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// End namespace Foam

// ************************************************************************* //
