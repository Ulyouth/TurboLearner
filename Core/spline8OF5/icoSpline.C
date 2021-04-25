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
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//namespace Foam
//{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
void Foam::icoSpline<Specie>::checkInputData() const
{
    if (Tlow_ >= Thigh_)
    {
        FatalErrorIn("icoSpline<Specie>::check()")
            << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }

    if (T1_ <= Tlow_)
    {
        FatalErrorIn("icoSpline<Specie>::check()")
            << "T1(" << T1_ << ") <= Tlow(" << Tlow_ << ')'
            << exit(FatalError);
    }

    if (T16_ > Thigh_)
    {
        FatalErrorIn("icoSpline<Specie>::check()")
            << "T16(" << T16_ << ") > Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::icoSpline<Specie>::icoSpline(Istream& is)
:
    Specie(is),
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
    T15_(readScalar(is)),
    T16_(readScalar(is))
    
{
    checkInputData();

    forAll(arhoCoeffs_, i)
    {
        is >> arhoCoeffs_[i];
    }

    forAll(brhoCoeffs_, i)
    {
        is >> brhoCoeffs_[i];
    }

    forAll(crhoCoeffs_, i)
    {
        is >> crhoCoeffs_[i];
    }

    forAll(drhoCoeffs_, i)
    {
        is >> drhoCoeffs_[i];
    }
    forAll(erhoCoeffs_, i)
    {
        is >> erhoCoeffs_[i];
    }
    forAll(frhoCoeffs_, i)
    {
        is >> frhoCoeffs_[i];
    }

    forAll(grhoCoeffs_, i)
    {
        is >> grhoCoeffs_[i];
    }

    forAll(hrhoCoeffs_, i)
    {
        is >> hrhoCoeffs_[i];
    }

    forAll(irhoCoeffs_, i)
    {
        is >> irhoCoeffs_[i];
    }
    forAll(jrhoCoeffs_, i)
    {
        is >> jrhoCoeffs_[i];
    }
    forAll(krhoCoeffs_, i)
    {
        is >> krhoCoeffs_[i];
    }

    forAll(lrhoCoeffs_, i)
    {
        is >> lrhoCoeffs_[i];
    }

    forAll(mrhoCoeffs_, i)
    {
        is >> mrhoCoeffs_[i];
    }

    forAll(nrhoCoeffs_, i)
    {
        is >> nrhoCoeffs_[i];
    }
    forAll(orhoCoeffs_, i)
    {
        is >> orhoCoeffs_[i];
    }
    forAll(prhoCoeffs_, i)
    {
        is >> prhoCoeffs_[i];
    }

    forAll(qrhoCoeffs_, i)
    {
        is >> qrhoCoeffs_[i];
    } 
    // Check state of Istream
    is.check("icoSpline::icoSpline(Istream& is)");
}

template<class Specie>
Foam::icoSpline<Specie>::icoSpline(const dictionary& dict)
:
    Specie(dict),
    Tlow_(readScalar(dict.subDict("equationOfState").lookup("Tlow"))),
    Thigh_(readScalar(dict.subDict("equationOfState").lookup("Thigh"))),
    T1_(readScalar(dict.subDict("equationOfState").lookup("T1"))),
    T2_(readScalar(dict.subDict("equationOfState").lookup("T2"))),
    T3_(readScalar(dict.subDict("equationOfState").lookup("T3"))),
    T4_(readScalar(dict.subDict("equationOfState").lookup("T4"))),
    T5_(readScalar(dict.subDict("equationOfState").lookup("T5"))),
    T6_(readScalar(dict.subDict("equationOfState").lookup("T6"))),
    T7_(readScalar(dict.subDict("equationOfState").lookup("T7"))),
    T8_(readScalar(dict.subDict("equationOfState").lookup("T8"))),
    T9_(readScalar(dict.subDict("equationOfState").lookup("T9"))),
    T10_(readScalar(dict.subDict("equationOfState").lookup("T10"))),
    T11_(readScalar(dict.subDict("equationOfState").lookup("T11"))),
    T12_(readScalar(dict.subDict("equationOfState").lookup("T12"))),
    T13_(readScalar(dict.subDict("equationOfState").lookup("T13"))),
    T14_(readScalar(dict.subDict("equationOfState").lookup("T14"))),
    T15_(readScalar(dict.subDict("equationOfState").lookup("T15"))),
    T16_(readScalar(dict.subDict("equationOfState").lookup("T16"))),
    arhoCoeffs_(dict.subDict("equationOfState").lookup("arhoCoeffs")),
    brhoCoeffs_(dict.subDict("equationOfState").lookup("brhoCoeffs")),
    crhoCoeffs_(dict.subDict("equationOfState").lookup("crhoCoeffs")),
    drhoCoeffs_(dict.subDict("equationOfState").lookup("drhoCoeffs")),
    erhoCoeffs_(dict.subDict("equationOfState").lookup("erhoCoeffs")),
    frhoCoeffs_(dict.subDict("equationOfState").lookup("frhoCoeffs")),
    grhoCoeffs_(dict.subDict("equationOfState").lookup("grhoCoeffs")),
    hrhoCoeffs_(dict.subDict("equationOfState").lookup("hrhoCoeffs")),
    irhoCoeffs_(dict.subDict("equationOfState").lookup("irhoCoeffs")),
    jrhoCoeffs_(dict.subDict("equationOfState").lookup("jrhoCoeffs")),
    krhoCoeffs_(dict.subDict("equationOfState").lookup("krhoCoeffs")),
    lrhoCoeffs_(dict.subDict("equationOfState").lookup("lrhoCoeffs")),
    mrhoCoeffs_(dict.subDict("equationOfState").lookup("mrhoCoeffs")),
    nrhoCoeffs_(dict.subDict("equationOfState").lookup("nrhoCoeffs")),
    orhoCoeffs_(dict.subDict("equationOfState").lookup("orhoCoeffs")),
    prhoCoeffs_(dict.subDict("equationOfState").lookup("prhoCoeffs")),
    qrhoCoeffs_(dict.subDict("equationOfState").lookup("qrhoCoeffs"))

{
    checkInputData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::icoSpline<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add("Tlow", Tlow_);
    dict.add("Thigh", Thigh_);
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
    dict.add("arhoCoeffs", arhoCoeffs_);
    dict.add("brhoCoeffs", brhoCoeffs_);
    dict.add("crhoCoeffs", crhoCoeffs_);
    dict.add("drhoCoeffs", drhoCoeffs_);
    dict.add("erhoCoeffs", erhoCoeffs_);
    dict.add("frhoCoeffs", frhoCoeffs_);
    dict.add("grhoCoeffs", grhoCoeffs_);
    dict.add("hrhoCoeffs", hrhoCoeffs_);
    dict.add("irhoCoeffs", irhoCoeffs_);
    dict.add("jrhoCoeffs", jrhoCoeffs_);
    dict.add("krhoCoeffs", krhoCoeffs_);
    dict.add("lrhoCoeffs", lrhoCoeffs_);
    dict.add("mrhoCoeffs", mrhoCoeffs_);
    dict.add("nrhoCoeffs", nrhoCoeffs_);
    dict.add("orhoCoeffs", orhoCoeffs_);
    dict.add("prhoCoeffs", prhoCoeffs_);
    dict.add("qrhoCoeffs", qrhoCoeffs_);

    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const icoSpline<Specie>& iS)
        
{
    os  << static_cast<const Specie&>(iS) << nl
        << "    " << iS.Tlow_
        << tab << iS.Thigh_
        << tab << iS.T1_
        << tab << iS.T2_
        << tab << iS.T3_
        << tab << iS.T4_
        << tab << iS.T5_
        << tab << iS.T6_
        << tab << iS.T7_
        << tab << iS.T8_
        << tab << iS.T9_
        << tab << iS.T10_
        << tab << iS.T11_
        << tab << iS.T12_
        << tab << iS.T13_
        << tab << iS.T14_
        << tab << iS.T15_
        << tab << iS.T16_;

    os << nl << "    ";

    forAll(iS.arhoCoeffs_, i)
    {
        os << iS.arhoCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.brhoCoeffs_, i)
    {
        os << iS.brhoCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.crhoCoeffs_, i)
    {
        os << iS.crhoCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.drhoCoeffs_, i)
    {
        os << iS.drhoCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.erhoCoeffs_, i)
    {
        os << iS.erhoCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.frhoCoeffs_, i)
    {
        os << iS.frhoCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.grhoCoeffs_, i)
    {
        os << iS.grhoCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.hrhoCoeffs_, i)
    {
        os << iS.hrhoCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.irhoCoeffs_, i)
    {
        os << iS.irhoCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.jrhoCoeffs_, i)
    {
        os << iS.jrhoCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.krhoCoeffs_, i)
    {
        os << iS.krhoCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.lrhoCoeffs_, i)
    {
        os << iS.lrhoCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.mrhoCoeffs_, i)
    {
        os << iS.mrhoCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.nrhoCoeffs_, i)
    {
        os << iS.nrhoCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.orhoCoeffs_, i)
    {
        os << iS.orhoCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.prhoCoeffs_, i)
    {
        os << iS.prhoCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(iS.qrhoCoeffs_, i)
    {
        os << iS.qrhoCoeffs_[i] << ' ';
    }


    os.check
    (
        "Ostream& operator<<(Ostream& os, const icoSpline<Specie>& iS)"
        //"Ostream& operator<<"
        //"(Ostream& os, const icoSpline<Specie>& iS)"
    );

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// End namespace Foam

// ************************************************************************* //
