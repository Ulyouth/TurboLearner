/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) w
2011-2016 OpenFOAM Foundation
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
#include "myMyongKasagiKE.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"    		//for y calculation
#include "wallFvPatch.H" 		//for fvWallPatch to find wall patches
#include "error.H"
#include "uniformDimensionedFields.H"
#include "dimensionedSymmTensor.H" 	//validate Gk
#include "dimensionedTensor.H" 		//validate Gk
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
namespace RASModels
{
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void myMyongKasagiKE<BasicTurbulenceModel>::myYPlus()	//Berechnung des y+ Wertes
{


	//WallShearStress
        const volVectorField& C = this->mesh_.C(); 				//Cell center: including internalField and boundaryField; actually slicedVolVectorField
        nu_ = this->nu();
	Reff_ = this->nu()*this->rho()*(dev(twoSymm(fvc::grad(this->U_))));
	Rxz_ = Reff_.component(symmTensor::XZ);

	const volVectorField& U = this->U_;					//Geschwindigkeiten ausschreiben
	const volScalarField Uz = U.component(vector::Z);			// z Komponente der Geschwindigkeit

	labelList wallPatchIDs;							// ?????????????????

	//iterate through wall patches
	forAll(C.boundaryField(),patchi)
	{
    	if (isA<wallFvPatch>(this->mesh_.boundary()[patchi]))
    		{
        	label wallPatchi = this->mesh_.boundaryMesh().findPatchID(this->mesh_.boundary()[patchi].name());
		wallPatchIDs.append(wallPatchi);
     		}
	}
	label wallCellID = 0;
	label wallPatchID = wallPatchIDs[wallCellID];
	label wallPatchIDsSize = wallPatchIDs.size();

//iterate through internal field of C
        forAll(C.internalField(),cellID) //cellID from 0 to 5099 //parallel -n4 0-1274
	{
                //find correct wallPatchID and wallCellID for rhoWall, nuWall and tauWall
	        label foundMatchingZ = 0;
		forAll(wallPatchIDs,i)
		{
		    scalarField zValuesOfWallPatch = C.boundaryField()[wallPatchIDs[i]].component(vector::Z);
		    if((zValuesOfWallPatch.size() == 0) && (wallPatchID != wallPatchIDs[wallPatchIDsSize-1]))		//???????????????????
		    //works only for two wall patches well. It adds up every loop till it reaches the last wall patch
		    {
			 wallPatchID = wallPatchIDs[i+1];
		    }
		    if((zValuesOfWallPatch.size() != 0) && (wallPatchID == wallPatchIDs[i]))
			//first term: when mpirun there are some processors with empy wall patch z value list
			//second term: search in one patch and if nothing found, don't search it again
		    {
		        if((C[cellID].z() < (zValuesOfWallPatch[wallCellID]+2*SMALL)) &&
			   (C[cellID].z() > (zValuesOfWallPatch[wallCellID]-2*SMALL)))
		        {
			    wallPatchID = wallPatchIDs[i];
			    foundMatchingZ = 1;
			    break;
		        }
			else if(wallCellID < (zValuesOfWallPatch.size()-1) &&
			       (C[cellID].z() < (zValuesOfWallPatch[wallCellID+1]+2*SMALL)) &&
                               (C[cellID].z() > (zValuesOfWallPatch[wallCellID+1]-2*SMALL))) //here could a segmentation fault happen, when wallCellID+1 exceeds the range and is still evaluated
			{
                            wallCellID = wallCellID+1;
                            wallPatchID = wallPatchIDs[i];
                            foundMatchingZ = 1;
			    break;
		        }
		    	else if((foundMatchingZ == 0) &&
				(wallCellID < (zValuesOfWallPatch.size()-1)))
		    	{

			    forAll(C.boundaryField()[wallPatchIDs[i]],wallCelli)
		    	    {
			        if((C[cellID].z() < (zValuesOfWallPatch[wallCelli]+2*SMALL)) &&
				   (C[cellID].z() > (zValuesOfWallPatch[wallCelli]-2*SMALL)))
			        {
			   	    wallCellID = wallCelli;
				    wallPatchID = wallPatchIDs[i];
				    foundMatchingZ = 1;
				    break;
			    	}
			    	else
			    	{
//				    Pout << "couldn't find " << C[cellID].z() << " in " << zValuesOfWallPatch << endl;
			    	}
		            }
		    	}
		      	else
		    	{
			    if(wallPatchID != wallPatchIDs[wallPatchIDsSize-1])
			    {
				wallPatchID = wallPatchIDs[i+1];
				wallCellID = 0;					// modifikation neu
			    }
		    	}
		    }
		}
                scalarField zValuesOfWallPatch = C.boundaryField()[wallPatchID].component(vector::Z);

                if((nu_.boundaryField()[wallPatchID][wallCellID] < ROOTVSMALL) &&
                   (nu_.boundaryField()[wallPatchID][wallCellID] > -ROOTVSMALL))
		{
		    Pout << "viscosity too small in yPlus calculation" << endl;
		    FatalErrorInFunction
        	        << "viscosity too small in yPlus calculation"
            	        << abort(FatalError);
		}
		else if((foundMatchingZ == 0) && (wallCellID < zValuesOfWallPatch.size()-1))
                {
//		    if(wallCellID == (zValuesOfWallPatch.size()-1)
		    Pout << "didn't find wall coordinate for yPlus calculation (cellID = " << cellID << endl;
		    Pout << "wallCellID = " << wallCellID << "; wallPatchID = " << wallPatchID << endl;
		    Pout << "zValuesOfWallPatch.size() = " << zValuesOfWallPatch.size() << endl;
		    Pout << "Proof: C[" << cellID << "].z(): " << C[cellID] << " == " << zValuesOfWallPatch[wallCellID] << endl;
                    printf("C[cellID].z() = %2.18f and zValuesOfWallPatch[wallCellID] = %2.18f \n",
			   C[cellID].z(),zValuesOfWallPatch[wallCellID]);

                    yPlus_[cellID] = -1;
                }
                else
		{
                    //calculate yPlus (semi local scaling) = sqrt(tauWall/rho)*y/nu
		   double Rxz = mag(Rxz_.boundaryField()[wallPatchID][wallCellID]);

	if(yPlus_Model_.value() == 1.0)
	{
                   yPlus_[cellID] = (::sqrt(Rxz/
                                           this->rho()[cellID])
                                    *y_[cellID])
                                    /nu_[cellID];
                   uPlus_[cellID] = Uz[cellID]/::sqrt(Rxz/
                                             this->rho()[cellID]);
	}
	else if(yPlus_Model_.value() == 2.0)
	{
		yPlus_[cellID] = 2.4*	sqrt(sqrt(k_[cellID])*y_[cellID]/nu_[cellID])
                                       +0.0025*pow(sqrt(k_[cellID])*y_[cellID]/nu_[cellID],2);
		uPlus_[cellID] = Uz[cellID]/::sqrt(Rxz/
                                             this->rho()[cellID]);
	}
	else if(yPlus_Model_.value() == 0.0)
	{
                   yPlus_[cellID] = (::sqrt(Rxz/
                                             this->rho().boundaryField()[wallPatchID][wallCellID])
                                      *y_[cellID])
                                      /nu_.boundaryField()[wallPatchID][wallCellID];

                   uPlus_[cellID] = Uz[cellID]/::sqrt(Rxz/
                                             this->rho().boundaryField()[wallPatchID][wallCellID]);
	}
	else
	{
	    Pout << "Choose a apropriate Model for yPlus" << endl;
            FatalErrorInFunction
            << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            << abort(FatalError);
	}
	}
	}

//iterate through boundaryField of C
        forAll(C.boundaryField(),patchi) //patchi from 0 to 6
        {
//	    Pout << "name of patch[" << patchi << "] = " << this->mesh_.boundary()[patchi].name() << endl;
	    if((this->mesh_.boundary()[patchi].type() == "wedge") ||
	       (this->mesh_.boundary()[patchi].type() == "empty"))
	    {
	    }
	    else
	    {
		if (isA<wallFvPatch>(this->mesh_.boundary()[patchi]))
		{
//		    Pout << "patch = " << this->mesh_.boundary()[patchi].name() << endl;
		    forAll(C.boundaryField()[patchi],cellID) //depending on patch 0-5099 or 0-1 or 0-100
		    //if a rank has zero entries in tankwall1, it just skips this loop
		    {
			if((nu_.boundaryField()[patchi][cellID] < ROOTVSMALL) &&
			   (nu_.boundaryField()[patchi][cellID] > -ROOTVSMALL))
                	{
                    	    Pout << "viscosity too small in yPlus calculation" << endl;
                    	    FatalErrorInFunction
                        	<< "viscosity too small in yPlus calculation"
                        	<< abort(FatalError);
                	}
                	else
                	{
                            //calculate yPlus (semi local scaling) for wall boundary
                  double Rxz = mag(Rxz_.boundaryField()[patchi][cellID]);
//                   if(Rxz < 0.03)
//                   {
//                        Rxz = 0.03;
//+++++                        Info << "tauWall corrected to 0.03" << endl;
//                   }
//                   else if(Rxz > 3)
//                   {
//                        Rxz = 3;
//                        Info << "tauWall corrected to 3" << endl;
//                   }
		if(yPlus_Model_.value() == 1.0)
		{

                            yPlus_.boundaryFieldRef()[patchi][cellID] =
				(::sqrt(Rxz/
                             this->rho().boundaryField()[patchi][cellID])
                            	*y_.boundaryField()[patchi][cellID])
                            	/nu_.boundaryField()[patchi][cellID];
		}
		else if(yPlus_Model_.value() == 2.0)
		{
                        yPlus_.boundaryFieldRef()[patchi][cellID] = 
					 2.4*sqrt(sqrt(k_.boundaryField()[patchi][cellID])*y_.boundaryField()[patchi][cellID]/nu_.boundaryField()[patchi][cellID])
                                       +0.0025*pow(sqrt(k_.boundaryField()[patchi][cellID])*y_.boundaryField()[patchi][cellID]/nu_.boundaryField()[patchi][cellID],2);
		}
		else if(yPlus_Model_.value() == 0.0)
		{
                            yPlus_.boundaryFieldRef()[patchi][cellID]=
				(::sqrt(Rxz/
                                 this->rho().boundaryField()[patchi][cellID])
                                *y_.boundaryField()[patchi][cellID])
                                /nu_.boundaryField()[patchi][cellID];
		}
	        else
        	{
        	Pout << "Choose a apropriate Model for yPlus" << endl;
            	FatalErrorInFunction
            	<< "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            	<< abort(FatalError);
        	}
		}
		}
	    	}
		else
		{
	            label wallCellID = 0;
        	    label wallPatchID = wallPatchIDs[wallCellID];

                    forAll(C.boundaryField()[patchi],cellID) //depending on patch 0-5099 or 0-1 or 0-100
                    {
			if(this->mesh_.boundary()[patchi].name() == "inlet")
			{
			    //take closest wall cell to inlet (here: z=0.0057125)
			    //take first wall (wallPatchIDs[0]) and first cell in it (0)
			    wallPatchID = wallPatchIDs[0]; //careful. This is only true for the geometry in this study
			    wallCellID = 0; //careful. This is only true for the geometry in this study (inlet is located at the bottom)

			    if((nu_.boundaryField()[patchi][cellID] < ROOTVSMALL) &&
			       (nu_.boundaryField()[patchi][cellID] > -ROOTVSMALL))
			    {
				Pout << "viscosity too small in yPlus calculation" << endl;
				FatalErrorInFunction
				    << "viscosity too small in yPlus calculation"
				    << abort(FatalError);
			    }
			    else
			    {
                            	//calculate yPlus (semi local scaling) for non wall boundary
                  double Rxz = mag(Rxz_.boundaryField()[wallPatchID][wallCellID]);
//                   if(Rxz < 0.03)
//                   {
//                        Rxz = 0.03;
//                        Info << "tauWall corrected to 0.03" << endl;
//                   }
//                   else if(Rxz > 3)
//+++                   {
//                        Rxz = 3;
//                        Info << "wtauWall corrected to 3" << endl;
//                   }
		if(yPlus_Model_.value() == 1.0)
		{
                            yPlus_.boundaryFieldRef()[patchi][cellID] =
                                (::sqrt(Rxz/
                             this->rho().boundaryField()[patchi][cellID])
                                *y_.boundaryField()[patchi][cellID])
                                /nu_.boundaryField()[patchi][cellID];
		}
		else if(yPlus_Model_.value() == 2.0)
		{
                        yPlus_.boundaryFieldRef()[patchi][cellID] = 2.4*sqrt(sqrt(k_.boundaryField()[patchi][cellID])*y_.boundaryField()[patchi][cellID]/nu_.boundaryField()[patchi][cellID])
                                       +0.0025*pow(sqrt(k_.boundaryField()[patchi][cellID])*y_.boundaryField()[patchi][cellID]/nu_.boundaryField()[patchi][cellID],2);
		}
		else if(yPlus_Model_.value() == 0.0)
		{
                                yPlus_.boundaryFieldRef()[patchi][cellID] =
                                    (::sqrt(Rxz/
                                    this->rho().boundaryField()[wallPatchID][wallCellID])
                                    *y_.boundaryField()[patchi][cellID])
                                    /nu_.boundaryField()[wallPatchID][wallCellID];
		}
	    	else
       		{
            	Pout << "Choose a apropriate Model for yPlus" << endl;
            	FatalErrorInFunction
            	<< "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            	<< abort(FatalError);
        	}
		}
		}
			else if(this->mesh_.boundary()[patchi].name() == "outlet")
			{
			    //take closest wall cell to outlet (here: z=1.51535)
			    //take last wall (wallPatchIDs[wallPatchIDs.size()-1]) and
			    //last cell in it (C.boundaryField()[wallPatchIDs[wallPatchIDs.size()-1]].size()-1)
                            wallPatchID = wallPatchIDs[wallPatchIDsSize-1]; //careful. This is only true for the geometry in this study
                            wallCellID = C.boundaryField()[wallPatchIDs[wallPatchIDsSize-1]].size()-1; //careful.
				//This is only true for the geometry in this study (outlet is located at the top)

                            if((nu_.boundaryField()[patchi][cellID] < ROOTVSMALL) &&
                               (nu_.boundaryField()[patchi][cellID] > -ROOTVSMALL))
                            {
                                Pout << "viscosity too small in yPlus calculation" << endl;
                                FatalErrorInFunction
                                    << "viscosity too small in yPlus calculation"
                                    << abort(FatalError);
                            }
                            else
                            {
                           	//calculate yPlus (semi local scaling) for non wall boundary
                  double Rxz = mag(Rxz_.boundaryField()[wallPatchID][wallCellID]);
//                   if(Rxz < 0.03)
//                   {
//                        Rxz = 0.03;
//                        Info << "tauWall corrected to 0.03" << endl;
//                   }
//                   else if(Rxz > 3)
//                   {
//                        Rxz = 3;
//                        Info << "tauWall corrected to 3" << endl;
//                   }
		if(yPlus_Model_.value() == 1.0)
		{
                            yPlus_.boundaryFieldRef()[patchi][cellID] =
                                (::sqrt(Rxz/
                             this->rho().boundaryField()[patchi][cellID])
                                *y_.boundaryField()[patchi][cellID])
                                /nu_.boundaryField()[patchi][cellID];
		}

		else if(yPlus_Model_.value() == 2.0)
		{
                        yPlus_.boundaryFieldRef()[patchi][cellID] = 
				         2.4*sqrt(sqrt(k_.boundaryField()[patchi][cellID])*y_.boundaryField()[patchi][cellID]/nu_.boundaryField()[patchi][cellID])
                                       +0.0025*pow(sqrt(k_.boundaryField()[patchi][cellID])*y_.boundaryField()[patchi][cellID]/nu_.boundaryField()[patchi][cellID],2);
		}
		else if(yPlus_Model_.value() == 0.0)
		{
                                //calculate yPlus for non wall boundary
                                yPlus_.boundaryFieldRef()[patchi][cellID] =
                                    (::sqrt(Rxz/
                                    this->rho().boundaryField()[wallPatchID][wallCellID])
                                    *y_.boundaryField()[patchi][cellID])
                                    /nu_.boundaryField()[wallPatchID][wallCellID];
		}
	    	else
        	{
            	Pout << "Choose a apropriate Model for yPlus" << endl;
            	FatalErrorInFunction
            	<< "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" 
            	<< abort(FatalError);
        	}
		}
		}
//			else
//			{
//			    //code located in "yPlus_boundary_of_front_and_back.C". There is are only procBoundaryXtoY left.
//			    //if e.g. yPlus for front or back need to be calculated, it needs to be uncommented
//		    	}
		    }
	    	}
	    }
	}
	return;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> myMyongKasagiKE<BasicTurbulenceModel>::fMu() const
{
    volScalarField Ret = sqr(k_)/(this->nu()*epsilon_);
    forAll(Ret,i)
    {
	if((Ret[i] < ROOTVSMALL) &&
           (Ret[i] > -ROOTVSMALL))
    	{
	    Pout << "Ret[" << i << "]= " << Ret[i] << endl;
	    Pout << "turbulent Reynoldsnumber in fmu function too small" << endl;
	    FatalErrorInFunction
	    << "turbulent Reynoldsnumber in fmu function too small"
	    << abort(FatalError);
	}
    }
    return
	((scalar(1)-exp(-yPlus_/70))
	*min(2.0, (scalar(1)+3.45/sqrt(mag(sqr(k_)/(this->nu()*epsilon_))))));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> myMyongKasagiKE<BasicTurbulenceModel>::f2() const
{
    return
	((scalar(1)-(2.0/9.0)*exp(-sqr((sqr(k_)/(this->nu()*epsilon_))/6)))
	*sqr(scalar(1)-exp(-yPlus_/5)));
}

template<class BasicTurbulenceModel>                                             // new COde +++++++++++++++++$
void myMyongKasagiKE<BasicTurbulenceModel>::f1_Prt()
{
        tmp<volTensorField> tgradU = fvc::grad(this->U_);
        const volScalarField UgradXZ = tgradU().component(symmTensor::XZ);	// Komponente du/dx
 
	const volScalarField myT = this->mesh_.objectRegistry::template
                                lookupObject<volScalarField>("T");
        tmp<volVectorField> tgradT = fvc::grad(myT);				//Vektor dt/dx
        const volScalarField TgradX = tgradT().component(vector::X);    	// Komponenten dt/dx

        const volScalarField myCp = this->mesh_.objectRegistry::template	// grad Cp
                                lookupObject<volScalarField>("Cp");
        tmp<volVectorField> tgradCp = fvc::grad(myCp);
	const volScalarField CpgradX = tgradCp().component(vector::X);

        const volScalarField myRho = this->mesh_.objectRegistry::template        // grad Rho
                                lookupObject<volScalarField>("rho");
        tmp<volVectorField> tgradRho = fvc::grad(myRho);
        const volScalarField RhogradX = tgradRho().component(vector::X);

        const uniformDimensionedVectorField& myg_ = this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

        const volScalarField myAlpha = this->mesh_.objectRegistry::template        //   alpha
                                lookupObject<volScalarField>("thermo:alpha");

	const volScalarField myNu = this->nu();        //   nu
	const volScalarField myNut = this->nut_;
	const volVectorField& C = this->mesh_.C();
	
//++++++++++++++++++++++++++++++++++++++++++ konstant ++++++++++++++++++++++++++++++++++++++++++++
	if(Prt_Model_.value() == 0.0)
	{ 
         myPrt_ = max(min(Uz_grad_min_const_/mag(fvc::grad(this->U_)), 4.0),4.0);
	}
//++++++++++++++++++++++++++++++++++++++++++ nach Bae +++++++++++++++++++++++++++++++++++++++++++
        else if(Prt_Model_.value() == 1.0)
        {
        myPrt_= 0.85-0.85*(scalar(1)-(1.0+min(mag(this->U_)/myRho*mag(RhogradX/max(mag(UgradXZ),Uz_grad_min_const_)),1.0))/
                (1.0+myT/myRho*mag(RhogradX/max(mag(TgradX),T_grad_min_const_))+myT/myCp*mag(CpgradX/max(mag(TgradX),T_grad_min_const_))))*
                0.5*(scalar(1)+tanh((50-yPlus_)/10))*(scalar(1)-exp(-yPlus_/70));
        }
//++++++++++++++++++++++++++++++++++++++++++ nach Oost +++++++++++++++++++++++++++++++++++++++++++
        else if(Prt_Model_.value() == 2.0)
	{
	 myPrt_ = min(0.85*(scalar(1)+beta_*mag(myg_)*(mag(TgradX))/(max(mag(UgradXZ*UgradXZ),100*Uz_grad_min_const_*Uz_grad_min_const_))),10.0) ;
	}
//++++++++++++++++++++++++++++++++++++++++++ nach Tien +++++++++++++++++++++++++++++++++++++++++++
	else if(Prt_Model_.value() == 3.0)
	{
	myPrt_ = max(min(Uz_grad_min_const_/mag(fvc::grad(this->U_)), 0.85),0.85);

	forAll(C.internalField(),j)
	{
		if (yPlus_[j] <= 5.0)
 		{
		myPrt_[j]=1.0;
		}
		else if ((yPlus_[j]>5.0)&&(yPlus_[j]<=100.0))
                {
                myPrt_[j]=scalar(0.5)+scalar(0.1)*(myNu[j]*myRho[j]/myAlpha[j])/max((myNut[j]/myNu[j]),0.001);
                }
		else if	( yPlus_[j] > 100.0)
		{
		myPrt_[j]=0.85;
		}
	}}
//++++++++++++++++++++++++++++++++++++++++++ nach TWL +++++++++++++++++++++++++++++++++++++++++++
	else if(Prt_Model_.value() == 4.0)
	{
	myPrt_ = max(min(Uz_grad_min_const_/mag(fvc::grad(this->U_)), 0.85),0.85);

        forAll(C.internalField(),j)
	{
                if ( myNut[j]/myNu[j] <= 0.2)
                {
                myPrt_[j]=1.0;
                }
                else if ((myNut[j]/myNu[j]>0.2)&&(myNut[j]/myNu[j]<=10.0))
                {
                myPrt_[j]=scalar(0.85)+(myNu[j]*myRho[j]/myAlpha[j])/10.0;
                }
                else if ( myNut[j]/myNu[j] > 10)
                {
                myPrt_[j]=0.85;
                }
	}}

	else
	{
            Pout << "Choose a apropriate Model for trubulent Prandtl number" << endl;
            FatalErrorInFunction
            << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            << abort(FatalError);
        }

return;
}

template<class BasicTurbulenceModel>
void myMyongKasagiKE<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*fMu()*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> myMyongKasagiKE<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> myMyongKasagiKE<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
myMyongKasagiKE<BasicTurbulenceModel>::myMyongKasagiKE
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    Prt_Model_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt_Model",
            this->coeffDict_,
            0
        )
    ),
    yPlus_Model_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "yPlus_Model",
            this->coeffDict_,
            0
        )
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.4
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.80
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            -0.33
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.4
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    Reff_
    (
        IOobject
        (
            "Reff",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	this->devRhoReff()
     ),

    Rxz_
    (
        IOobject
        (
            "Rxz",
             this->runTime_.timeName(),
             this->mesh_,
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
	 this->mesh_,
         dimensionedScalar("Rxz",dimMass/(dimLength*sqr(dimTime)), 0.0)
    ),

    nu_
    (
	IOobject
	(
            "nu",
            this->runTime_.timeName(),
            this->mesh_,
	    IOobject::NO_READ,
            IOobject::NO_WRITE
	),
        this->mesh_,
        dimensionedScalar("nu",sqr(dimLength)/dimTime, 0.0)
    ),

    myNut_
    (
	IOobject
	(
            "myNut",
            this->runTime_.timeName(),
            this->mesh_,
	    IOobject::NO_READ,
            IOobject::NO_WRITE
	),
        this->mesh_,
        dimensionedScalar("myNut",sqr(dimLength)/dimTime, 0.0)
    ),

    uPlus_
    (
        IOobject
        (
            "uPlus",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("no dimension", dimless, 0.0)
    ),

    Ret_
    (
        IOobject
        (
            "Ret",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("no dimension", dimless, 0.0)
    ),

    beta_
    (
        IOobject
        (
            "beta",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("no dimension", dimMass/(dimMass*dimTemperature), 0.0)
    ),

    T_grad_min_const_		// new
    (				// <- Diese Variable wird verwendet, um zu vermeiden, dass bei der Pr_t Berechnung der Nenner nicht
                                // <- null wird und dass die Einheiten in der min-max-Funktion stimmen. Somit wird durch einen sehr kleinen Wert geteilt.
        IOobject
        (
            "T_grad_min_const_",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("no dimension", dimensionSet(0, -1, 0, 1, 0, 0, 0) , 0.0001)
    ),

    Uz_grad_min_const_         // new
    (			       // <- Diese Variable wird verwendet, um zu vermeiden, dass bei der Pr_t Berechnung der Nenner nicht
                               // <- null wird und dass die Einheiten in der min-max-Funktion stimmen. Somit wird durch einen sehr kleinen Wert geteilt.
        IOobject
        (
            "Uz_grad_min_const",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("no dimension", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0001)
    ),

    product_OF_eps_
    (
        IOobject
        (
            "product_OF_eps",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("product_OF_eps", dimMass/(dimLength*sqr(dimTime)*sqr(dimTime)), 0.0)
    ),

    product_ass_eps_
    (
        IOobject
        (
            "product_ass_eps",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("product_ass_eps", dimMass/(dimLength*sqr(dimTime)*sqr(dimTime)), 0.0)
    ),

    buoyancy_OF_eps_
    (
        IOobject
        (
            "buoyancy_OF_eps",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("buoyancy_OF_eps", dimMass/(dimLength*sqr(dimTime)*sqr(dimTime)), 0.0)
    ),

    buoyancy_ass_eps_
    (
        IOobject
        (
            "buoyancy_ass_eps",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("buoyancy_ass_eps", dimMass/(dimLength*sqr(dimTime)*sqr(dimTime)), 0.0)
    ),

    product_OF_k_
    (
        IOobject
        (
            "product_OF_k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("product_OF_k", dimMass/(dimLength*dimTime*sqr(dimTime)), 0.0)
    ),

    product_ass_k_
    (
        IOobject
        (
            "product_ass_k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("product_ass_k", dimMass/(dimLength*dimTime*sqr(dimTime)), 0.0)
    ),

    buoyancy_OF_k_
    (
        IOobject
        (
            "buoyancy_OF_k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("buoyancy_OF_k", dimMass/(dimLength*dimTime*sqr(dimTime)), 0.0)
    ),

    buoyancy_ass_k_
    (
        IOobject
        (
            "buoyancy_ass_k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("buoyancy_ass_k", dimMass/(dimLength*dimTime*sqr(dimTime)), 0.0)
    ),

    dissip_k_
    (
        IOobject
        (
            "dissip_k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("dissip_k", dimMass/(dimLength*dimTime*sqr(dimTime)), 0.0)
    ),

    dissip_eps_
    (
        IOobject
        (
            "dissip_eps",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("dissip_eps", dimMass/(dimLength*sqr(dimTime)*sqr(dimTime)), 0.0)
    ),

    myfMu_
    (
        IOobject
        (
            "myfMu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("no dimension", dimless, 0.0)
    ),

    myrho_
    (
        IOobject
        (
            "myrho",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("myrho", dimMass/(dimLength*sqr(dimLength)), 0.0)
    ),

    myf2_
    (
        IOobject
        (
            "myf2",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("no dimension", dimless, 0.0)
    ),


//   myf1_Prt_                                      // new COde +++++++++++++++++++++++++++++++++++++++++
//(
//        IOobject
//        (
//            "myf1_Prt",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        this->mesh_,
//        dimensionedScalar("no dimension", dimless, 0.0)
//    ),

    myPrt_                                      // new COde ++++++++++++++++++++++++++++++++++++++++++++++++++$
    (
        IOobject
        (
            "myPrt",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("no dimension", dimless, 0.0)
    ),



    yPlus_
    (
        IOobject
        (
            "yPlus",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("no dimension", dimless, 0.0)
    ),

    y_(wallDist::New(this->mesh_).y())
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool myMyongKasagiKE<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void myMyongKasagiKE<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    //calculate yPlus
    myYPlus();

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    myrho_ = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& myT = this->mesh_.objectRegistry::template
			        lookupObject<volScalarField>("T");
    const volVectorField& C = this->mesh_.C();
    volScalarField& nut = this->nut_;
    myNut_ = nut;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    // Calculate parameters and coefficients for Myong Kasagi low-Reynolds
    // number model

    Ret_ = sqr(k_)/(this->nu()*epsilon_);

    tmp<volTensorField> tgradU = fvc::grad(U);
    tmp<volVectorField> tgradT = fvc::grad(myT);

    myfMu_ = fMu();
    myf2_ = f2();
//    myf1_Prt_ = f1_Prt();

    volScalarField G(this->GName(), nut*(tgradU() && dev(twoSymm(tgradU()))));

// calculate beta as function of temperature
	forAll(C.internalField(),cellID)
	{

if(myT[cellID]<=220.15)
{
   beta_[cellID]=0.0025; 
 }
 else if((myT[cellID]>220.15)&&(myT[cellID]<=263.15)) 
{ 
 beta_[cellID]=0.0000000000e+00*pow((myT[cellID]-273.15),3)+1.7513486438e-07*pow((myT[cellID]-273.15),2)+5.1203263735e-05*pow((myT[cellID]-273.15),1) + 5.1318191509e-03;
}
 else if((myT[cellID]>263.15)&&(myT[cellID]<=283.15)) 
{ 
 beta_[cellID]=1.9284184076e-07*pow((myT[cellID]-273.15),3)+5.9603900843e-06*pow((myT[cellID]-273.15),2)+1.0905581593e-04*pow((myT[cellID]-273.15),1) + 5.3246609916e-03;
}
 else if((myT[cellID]>283.15)&&(myT[cellID]<=293.15)) 
{ 
 beta_[cellID]=-4.1382622886e-07*pow((myT[cellID]-273.15),3)+2.4160432173e-05*pow((myT[cellID]-273.15),2)+-7.2944604942e-05*pow((myT[cellID]-273.15),1) + 5.9313290611e-03;
}
 else if((myT[cellID]>293.15)&&(myT[cellID]<=303.15)) 
{ 
 beta_[cellID]=1.3825102083e-05*pow((myT[cellID]-273.15),3)+-8.3017526653e-04*pow((myT[cellID]-273.15),2)+1.7013769369e-02*pow((myT[cellID]-273.15),1) + -1.0798009743e-01;
}
 else if((myT[cellID]>303.15)&&(myT[cellID]<=305.15)) 
{ 
 beta_[cellID]=1.1647007893e-03*pow((myT[cellID]-273.15),3)+-1.0440898712e-01*pow((myT[cellID]-273.15),2)+3.1243781248e+00*pow((myT[cellID]-273.15),1) + -3.1181623651e+01;
}
 else if((myT[cellID]>305.15)&&(myT[cellID]<=306.15)) 
{ 
 beta_[cellID]=1.3329371877e-03*pow((myT[cellID]-273.15),3)+-1.2055968135e-01*pow((myT[cellID]-273.15),2)+3.6412003398e+00*pow((myT[cellID]-273.15),1) + -3.6694393947e+01;
}
 else if((myT[cellID]>306.15)&&(myT[cellID]<=307.15)) 
{ 
 beta_[cellID]=5.3243736749e-02*pow((myT[cellID]-273.15),3)+-5.2597288380e+00*pow((myT[cellID]-273.15),2)+1.7323378251e+02*pow((myT[cellID]-273.15),1) + -1.9022127978e+03;
}
 else if((myT[cellID]>307.15)&&(myT[cellID]<=307.85)) 
{ 
 beta_[cellID]=-3.5208522953e-01*pow((myT[cellID]-273.15),3)+3.6083825722e+01*pow((myT[cellID]-273.15),2)+-1.2324470725e+03*pow((myT[cellID]-273.15),1) + 1.4028836893e+04;
}
 else if((myT[cellID]>307.85)&&(myT[cellID]<=308.15)) 
{ 
 beta_[cellID]=7.1503867240e-01*pow((myT[cellID]-273.15),3)+-7.5003772469e+01*pow((myT[cellID]-273.15),2)+2.6222925847e+03*pow((myT[cellID]-273.15),1) + -3.0557651809e+04;
}
 else if((myT[cellID]>308.15)&&(myT[cellID]<=309.15)) 
{ 
 beta_[cellID]=-7.5010931995e-03*pow((myT[cellID]-273.15),3)+8.6290291956e-01*pow((myT[cellID]-273.15),2)+-3.3041053902e+01*pow((myT[cellID]-273.15),1) + 4.2124064102e+02;
}
 else if((myT[cellID]>309.15)&&(myT[cellID]<=310.15)) 
{ 
 beta_[cellID]=-1.6867801235e-02*pow((myT[cellID]-273.15),3)+1.8745073874e+00*pow((myT[cellID]-273.15),2)+-6.9458814746e+01*pow((myT[cellID]-273.15),1) + 8.5825377115e+02;
}
 else if((myT[cellID]>310.15)&&(myT[cellID]<=311.15)) 
{ 
 beta_[cellID]=7.1329814131e-04*pow((myT[cellID]-273.15),3)+-7.6994643373e-02*pow((myT[cellID]-273.15),2)+2.7467603948e+00*pow((myT[cellID]-273.15),1) + -3.2281655582e+01;
}
 else if((myT[cellID]>311.15)&&(myT[cellID]<=313.15)) 
{ 
 beta_[cellID]=-6.3030901629e-04*pow((myT[cellID]-273.15),3)+7.6176572593e-02*pow((myT[cellID]-273.15),2)+-3.0737458119e+00*pow((myT[cellID]-273.15),1) + 4.1444756369e+01;
}
 else if((myT[cellID]>313.15)&&(myT[cellID]<=318.15)) 
{ 
 beta_[cellID]=-2.6630821176e-05*pow((myT[cellID]-273.15),3)+3.7351891789e-03*pow((myT[cellID]-273.15),2)+-1.7609047533e-01*pow((myT[cellID]-273.15),1) + 2.8093518821e+00;
}
 else if((myT[cellID]>318.15)&&(myT[cellID]<=323.15)) 
{ 
 beta_[cellID]=-6.7541492180e-06*pow((myT[cellID]-273.15),3)+1.0518384645e-03*pow((myT[cellID]-273.15),2)+-5.5339693179e-02*pow((myT[cellID]-273.15),1) + 9.9809014986e-01;
}
 else if((myT[cellID]>323.15)&&(myT[cellID]<=333.15)) 
{ 
 beta_[cellID]=-9.6294957959e-07*pow((myT[cellID]-273.15),3)+1.8315851877e-04*pow((myT[cellID]-273.15),2)+-1.1905695891e-02*pow((myT[cellID]-273.15),1) + 2.7419019505e-01;
}
 else if((myT[cellID]>333.15)&&(myT[cellID]<=343.15)) 
{ 
 beta_[cellID]=-1.8786846884e-07*pow((myT[cellID]-273.15),3)+4.3643918838e-05*pow((myT[cellID]-273.15),2)+-3.5348198946e-03*pow((myT[cellID]-273.15),1) + 1.0677267513e-01;
}
 else if((myT[cellID]>343.15)&&(myT[cellID]<=363.15)) 
{ 
 beta_[cellID]=-5.2220145788e-08*pow((myT[cellID]-273.15),3)+1.5157770996e-05*pow((myT[cellID]-273.15),2)+-1.5407895457e-03*pow((myT[cellID]-273.15),1) + 6.0245300324e-02;
}
 else if((myT[cellID]>363.15)&&(myT[cellID]<=383.15)) 
{ 
 beta_[cellID]=-8.1408091640e-09*pow((myT[cellID]-273.15),3)+3.2563501080e-06*pow((myT[cellID]-273.15),2)+-4.6966166577e-04*pow((myT[cellID]-273.15),1) + 2.8111463925e-02;
}
 else if((myT[cellID]>383.15)&&(myT[cellID]<=403.15)) 
{ 
 beta_[cellID]=-4.8666175563e-09*pow((myT[cellID]-273.15),3)+2.1758668775e-06*pow((myT[cellID]-273.15),2)+-3.5080851041e-04*pow((myT[cellID]-273.15),1) + 2.3753514895e-02;
}
 else if((myT[cellID]>403.15)&&(myT[cellID]<=423.15)) 
{ 
 beta_[cellID]=-1.9802206109e-09*pow((myT[cellID]-273.15),3)+1.0501720688e-06*pow((myT[cellID]-273.15),2)+-2.0446818528e-04*pow((myT[cellID]-273.15),1) + 1.7412100806e-02;
}
 else if((myT[cellID]>423.15)&&(myT[cellID]<=443.15)) 
{ 
 beta_[cellID]=7.1523436044e-25*pow((myT[cellID]-273.15),3)+1.5907279389e-07*pow((myT[cellID]-273.15),2)+-7.0803294045e-05*pow((myT[cellID]-273.15),1) + 1.0728856244e-02;
}
else if(myT[cellID]>443.15)
{
   beta_[cellID]=0.002; 
}

	    else
	    {
                    Pout << "no matching beta for temperature in turbulence model" << endl;
                    FatalErrorInFunction
                        << "no matching beta for temperature in turbulence model"
                        << abort(FatalError);
	    }
	}

	forAll(C.boundaryField(),patchi) //patchi from 0 to 6
       {
	    forAll(C.boundaryField()[patchi],cellID)
            {

if(myT.boundaryField()[patchi][cellID]<=220.15)
{
   beta_.boundaryFieldRef()[patchi][cellID]=0.0025; 
 }
 else if((myT.boundaryField()[patchi][cellID]>220.15)&&(myT.boundaryField()[patchi][cellID]<=263.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=0.0000000000e+00*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+1.7513486438e-07*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+5.1203263735e-05*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+5.1318191509e-03;
}
 else if((myT.boundaryField()[patchi][cellID]>263.15)&&(myT.boundaryField()[patchi][cellID]<=283.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=1.9284184076e-07*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+5.9603900843e-06*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+1.0905581593e-04*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+5.3246609916e-03;
}
 else if((myT.boundaryField()[patchi][cellID]>283.15)&&(myT.boundaryField()[patchi][cellID]<=293.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=-4.1382622886e-07*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+2.4160432173e-05*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-7.2944604942e-05*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+5.9313290611e-03;
}
 else if((myT.boundaryField()[patchi][cellID]>293.15)&&(myT.boundaryField()[patchi][cellID]<=303.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=1.3825102083e-05*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+-8.3017526653e-04*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+1.7013769369e-02*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+-1.0798009743e-01;
}
 else if((myT.boundaryField()[patchi][cellID]>303.15)&&(myT.boundaryField()[patchi][cellID]<=305.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=1.1647007893e-03*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+-1.0440898712e-01*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+3.1243781248e+00*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+-3.1181623651e+01;
}
 else if((myT.boundaryField()[patchi][cellID]>305.15)&&(myT.boundaryField()[patchi][cellID]<=306.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=1.3329371877e-03*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+-1.2055968135e-01*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+3.6412003398e+00*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+-3.6694393947e+01;
}
 else if((myT.boundaryField()[patchi][cellID]>306.15)&&(myT.boundaryField()[patchi][cellID]<=307.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=5.3243736749e-02*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+-5.2597288380e+00*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+1.7323378251e+02*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+-1.9022127978e+03;
}
 else if((myT.boundaryField()[patchi][cellID]>307.15)&&(myT.boundaryField()[patchi][cellID]<=307.85)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=-3.5208522953e-01*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+3.6083825722e+01*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-1.2324470725e+03*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+1.4028836893e+04;
}
 else if((myT.boundaryField()[patchi][cellID]>307.85)&&(myT.boundaryField()[patchi][cellID]<=308.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=7.1503867240e-01*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+-7.5003772469e+01*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+2.6222925847e+03*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+-3.0557651809e+04;
}
 else if((myT.boundaryField()[patchi][cellID]>308.15)&&(myT.boundaryField()[patchi][cellID]<=309.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=-7.5010931995e-03*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+8.6290291956e-01*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-3.3041053902e+01*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+4.2124064102e+02;
}
 else if((myT.boundaryField()[patchi][cellID]>309.15)&&(myT.boundaryField()[patchi][cellID]<=310.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=-1.6867801235e-02*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+1.8745073874e+00*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-6.9458814746e+01*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+8.5825377115e+02;
}
 else if((myT.boundaryField()[patchi][cellID]>310.15)&&(myT.boundaryField()[patchi][cellID]<=311.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=7.1329814131e-04*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+-7.6994643373e-02*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+2.7467603948e+00*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+-3.2281655582e+01;
}
 else if((myT.boundaryField()[patchi][cellID]>311.15)&&(myT.boundaryField()[patchi][cellID]<=313.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=-6.3030901629e-04*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+7.6176572593e-02*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-3.0737458119e+00*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+4.1444756369e+01;
}
 else if((myT.boundaryField()[patchi][cellID]>313.15)&&(myT.boundaryField()[patchi][cellID]<=318.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=-2.6630821176e-05*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+3.7351891789e-03*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-1.7609047533e-01*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+2.8093518821e+00;
}
 else if((myT.boundaryField()[patchi][cellID]>318.15)&&(myT.boundaryField()[patchi][cellID]<=323.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=-6.7541492180e-06*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+1.0518384645e-03*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-5.5339693179e-02*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+9.9809014986e-01;
}
 else if((myT.boundaryField()[patchi][cellID]>323.15)&&(myT.boundaryField()[patchi][cellID]<=333.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=-9.6294957959e-07*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+1.8315851877e-04*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-1.1905695891e-02*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+2.7419019505e-01;
}
 else if((myT.boundaryField()[patchi][cellID]>333.15)&&(myT.boundaryField()[patchi][cellID]<=343.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=-1.8786846884e-07*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+4.3643918838e-05*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-3.5348198946e-03*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+1.0677267513e-01;
}
 else if((myT.boundaryField()[patchi][cellID]>343.15)&&(myT.boundaryField()[patchi][cellID]<=363.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=-5.2220145788e-08*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+1.5157770996e-05*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-1.5407895457e-03*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+6.0245300324e-02;
}
 else if((myT.boundaryField()[patchi][cellID]>363.15)&&(myT.boundaryField()[patchi][cellID]<=383.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=-8.1408091640e-09*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+3.2563501080e-06*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-4.6966166577e-04*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+2.8111463925e-02;
}
 else if((myT.boundaryField()[patchi][cellID]>383.15)&&(myT.boundaryField()[patchi][cellID]<=403.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=-4.8666175563e-09*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+2.1758668775e-06*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-3.5080851041e-04*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+2.3753514895e-02;
}
 else if((myT.boundaryField()[patchi][cellID]>403.15)&&(myT.boundaryField()[patchi][cellID]<=423.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=-1.9802206109e-09*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+1.0501720688e-06*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-2.0446818528e-04*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+1.7412100806e-02;
}
 else if((myT.boundaryField()[patchi][cellID]>423.15)&&(myT.boundaryField()[patchi][cellID]<=443.15)) 
{ 
 beta_.boundaryFieldRef()[patchi][cellID]=7.1523436044e-25*pow((myT.boundaryField()[patchi][cellID]-273.15),3)+1.5907279389e-07*pow((myT.boundaryField()[patchi][cellID]-273.15),2)+-7.0803294045e-05*pow((myT.boundaryField()[patchi][cellID]-273.15),1)+1.0728856244e-02;
}
else if(myT.boundaryField()[patchi][cellID]>443.15)
{
   beta_.boundaryFieldRef()[patchi][cellID]=0.002; 
}


		else
		{
                    Pout << "no matching beta for temperature in turbulence model" << endl;
                    FatalErrorInFunction
                        << "no matching beta for temperature in turbulence model"
                        << abort(FatalError);
		}
	    }
	}

        const uniformDimensionedVectorField& myg_ = this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");
        tensor I(1,0,0,0,1,0,0,0,1);
	volScalarField P_k = nut*(tgradU() && twoSymm(tgradU()));    // dev !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	volScalarField G_k = (-0.3)*beta_*k_/epsilon_*
			     (myg_ & ((nut*rho*twoSymm(tgradU())-2.0/3.0*rho*k_*I) & tgradT()));
	scalar minus = -1.0;

	product_OF_eps_ = C1_*alpha*rho*G*epsilon_/k_;
	buoyancy_OF_eps_ = minus*((2.0/3.0)*C1_ + C3_)*alpha*rho*divU*epsilon_;

	product_ass_eps_ = rho*C1_*epsilon_/k_*P_k;
	buoyancy_ass_eps_ = C1_*epsilon_/k_*G_k;

	product_OF_k_ = alpha*rho*G;
	buoyancy_OF_k_ = minus*2.0/3.0*alpha*rho*divU*k_;

	product_ass_k_ = rho*P_k;
	buoyancy_ass_k_ = G_k;

	dissip_eps_ = minus*C2_*f2()*alpha*rho*epsilon_*epsilon_/k_;
	dissip_k_ = minus*alpha*rho*epsilon_;

//         myPrt_ =  f1_Prt();
	f1_Prt();
//        forAll(C.internalField(),cellID)
//        {
//	if((myT[cellID]>305.0)&&(myT[cellID]<=310.0))
//		{
//		myPrt_ = 3.0; 
//        	}
//	else
//		{
//		 myPrt_ = 0.85;
//		}
//	}
        //this->alphat()=this->mut()/myPrt_;

	tgradU.clear();
	tgradT.clear();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
	  rho*C1_*epsilon_/k_*P_k //production assessment (Assessment of low-Reynolds number k-epsilon models aigainst highly buoyant flows (Bae, Kim, Kim 2016)
//        - fvm::SuSp(((2.0/3.0)*C1_)*alpha*rho*divU, epsilon_) //new
        + fvm::SuSp(C1_/k_*G_k, epsilon_) //buoyancy assessment
        - fvm::Sp(C2_*f2()*alpha*rho*epsilon_/k_, epsilon_) // SuSp   dissipation assessment and OF
//      + epsilonSource()
//      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
	  rho*P_k //production assessment
//        - fvm::SuSp(2.0/3.0*alpha*rho*divU, k_)    // new 
	+ fvm::SuSp(G_k/k_, k_)	//buoyancy assessment
	- fvm::Sp(alpha*rho*(epsilon_)/k_, k_) //dissipation assessment and OF
//      + kSource()
//      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
    this->alphat_=this->mut()/myPrt_;     // Überschreibe von alphat
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
