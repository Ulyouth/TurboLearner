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

#include "heRhoThermonn.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermonn<BasicPsiThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he().internalField();
    const scalarField& pCells = this->p_.internalField();

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        if (hCells[celli]<174644.281885) //-10C
        { TCells[celli]=  -0.0000016455*pow((hCells[celli]/1000),3)+0.0003455352*pow((hCells[celli]/1000),2)+0.4975307317*pow((hCells[celli]/1000),1)-99.0025109166+273.15;}


        else if (hCells[celli]<220035.161007) //10C
        { TCells[celli]=  -0.0000033379*pow((hCells[celli]/1000),3)+0.0012361678*pow((hCells[celli]/1000),2)+0.3444104907*pow((hCells[celli]/1000),1)-90.0725385826+273.15;}


        else if (hCells[celli]<284035.445055) //30C
        { TCells[celli]=  -0.0000073652*pow((hCells[celli]/1000),3)+0.00392816*pow((hCells[celli]/1000),2)-0.256465731*pow((hCells[celli]/1000),1)-45.2876197102+273.15;}


        else if (hCells[celli]<296424.864074) //32C
        { TCells[celli]=  0.0000022256*pow((hCells[celli]/1000),3)-0.0043749198*pow((hCells[celli]/1000),2)+2.1384582394*pow((hCells[celli]/1000),1)-275.4450093257+273.15;}


        else if (hCells[celli]<320691.340900) //34C
        { TCells[celli]=  0.0000193091*pow((hCells[celli]/1000),3)-0.0196829216*pow((hCells[celli]/1000),2)+6.7110823703*pow((hCells[celli]/1000),1)-730.7671044576+273.15;}


        else if (hCells[celli]<352291.220422) //35C
        { TCells[celli]=  0.000013184675*pow((hCells[celli]/1000),3)-0.013438531271*pow((hCells[celli]/1000),2)+4.593640957965*pow((hCells[celli]/1000),1)-491.921898833989+273.15;}


        else if (hCells[celli]<383341.178844) //37C
        { TCells[celli]=  0.00001522049*pow((hCells[celli]/1000),3)-0.015579004264*pow((hCells[celli]/1000),2)+5.343689864658*pow((hCells[celli]/1000),1)-579.517774176279+273.15;}


        else if (hCells[celli]<402901.495719) //40C
        { TCells[celli]=  0.000013156113*pow((hCells[celli]/1000),3)-0.013143500293*pow((hCells[celli]/1000),2)+4.386493374288*pow((hCells[celli]/1000),1)-454.192156340538+273.15;}


        else if (hCells[celli]<506211.457844) //90C
        { TCells[celli]=  -0.000006021*pow((hCells[celli]/1000),3)+0.0107761544*pow((hCells[celli]/1000),2)-5.5648733469*pow((hCells[celli]/1000),1)+926.6637070459+273.15;}


        else if (hCells[celli]<605294.947744) //170C
        { TCells[celli]=  -0.0000038579*pow((hCells[celli]/1000),3)+0.0071717453*pow((hCells[celli]/1000),2)-3.5798833965*pow((hCells[celli]/1000),1)+564.8344923959+273.15;}


        else 

        { TCells[celli]=-0.0000007128*pow((hCells[celli]/1000),3)+0.0014932393*pow((hCells[celli]/1000),2)-0.1564421351*pow((hCells[celli]/1000),1)-124.3454641833+273.15;}


        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();	
	
    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& ph = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                ph[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);


        if (ph[facei]<174644.281885) //-10C
        { pT[facei]=  -0.0000016455*pow((ph[facei]/1000),3)+0.0003455352*pow((ph[facei]/1000),2)+0.4975307317*pow((ph[facei]/1000),1)-99.0025109166+273.15;}


        else if (ph[facei]<220035.161007) //10C
        { pT[facei]=  -0.0000033379*pow((ph[facei]/1000),3)+0.0012361678*pow((ph[facei]/1000),2)+0.3444104907*pow((ph[facei]/1000),1)-90.0725385826+273.15;}


        else if (ph[facei]<284035.445055) //30C
        { pT[facei]=  -0.0000073652*pow((ph[facei]/1000),3)+0.00392816*pow((ph[facei]/1000),2)-0.256465731*pow((ph[facei]/1000),1)-45.2876197102+273.15;}


        else if (ph[facei]<296424.864074) //32C
        { pT[facei]=  0.0000022256*pow((ph[facei]/1000),3)-0.0043749198*pow((ph[facei]/1000),2)+2.1384582394*pow((ph[facei]/1000),1)-275.4450093257+273.15;}


        else if (ph[facei]<320691.340900) //34C
        { pT[facei]=  0.0000193091*pow((ph[facei]/1000),3)-0.0196829216*pow((ph[facei]/1000),2)+6.7110823703*pow((ph[facei]/1000),1)-730.7671044576+273.15;}


        else if (ph[facei]<352291.220422) //35C
        { pT[facei]=  0.000013184675*pow((ph[facei]/1000),3)-0.013438531271*pow((ph[facei]/1000),2)+4.593640957965*pow((ph[facei]/1000),1)-491.921898833989+273.15;}


        else if (ph[facei]<383341.178844) //37C
        { pT[facei]=  0.00001522049*pow((ph[facei]/1000),3)-0.015579004264*pow((ph[facei]/1000),2)+5.343689864658*pow((ph[facei]/1000),1)-579.517774176279+273.15;}


        else if (ph[facei]<402901.495719) //40C
        { pT[facei]=  0.000013156113*pow((ph[facei]/1000),3)-0.013143500293*pow((ph[facei]/1000),2)+4.386493374288*pow((ph[facei]/1000),1)-454.192156340538+273.15;}


        else if (ph[facei]<506211.457844) //90C
        { pT[facei]=  -0.000006021*pow((ph[facei]/1000),3)+0.0107761544*pow((ph[facei]/1000),2)-5.5648733469*pow((ph[facei]/1000),1)+926.6637070459+273.15;}


        else if (ph[facei]<605294.947744) //170C
        { pT[facei]=  -0.0000038579*pow((ph[facei]/1000),3)+0.0071717453*pow((ph[facei]/1000),2)-3.5798833965*pow((ph[facei]/1000),1)+564.8344923959+273.15;}


        else 

        { pT[facei]=-0.0000007128*pow((ph[facei]/1000),3)+0.0014932393*pow((ph[facei]/1000),2)-0.1564421351*pow((ph[facei]/1000),1)-124.3454641833+273.15;}

                //pT[facei] = mixture_.THE(ph[facei], pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRhoThermonn<BasicPsiThermo, MixtureType>::heRhoThermonn
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRhoThermonn<BasicPsiThermo, MixtureType>::~heRhoThermonn()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermonn<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering heRhoThermonn<MixtureType>::correct()" << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "exiting heRhoThermonn<MixtureType>::correct()" << endl;
    }
}


// ************************************************************************* //
