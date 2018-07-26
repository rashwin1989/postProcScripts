/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    calcDropC10

Description
    Calculates total normalized C10 inside droplet at each time

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvc.H"
#include "OFstream.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"            

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fileName outputFile1("ist_vs_t");
    OFstream os(runTime.path()/outputFile1);

    //scalar pi = Foam::constant::mathematical::pi;

    //scalar dz;
    //dz = runTime.controlDict().lookupOrDefault("dz",0.0001);    

    scalar AInterface;
    scalar VolPhase0;
    scalar VolPhase1;
    scalar m00;
    scalar m01;
    scalar m10;
    scalar m11;   
    vector zeroVec(0,0,0);
    vector velWaterPhase;
    vector velOilPhase;
    scalar mDot00;
    scalar mDot01;
    scalar mDot10;
    scalar mDot11;


    const scalar rho0 = 275;
    const scalar rho1 = 900;
    const scalarField& V = mesh.V();    

    while (runTime.run())
    {        
        runTime++;        

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Reading fields

    volScalarField alpha1
        (
            IOobject
            (
                "alpha1",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        ); 

        volScalarField A_intfc
        (
            IOobject
            (
                "A_intfc",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        ); 

    volScalarField C10
        (
            IOobject
            (
                "C10",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        ); 


    volScalarField C00
        (
            IOobject
            (
                "C00",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        ); 


     volScalarField Y10
        (
            IOobject
            (
                "Y10",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        ); 

       
     volScalarField Y00
        (
            IOobject
            (
                "Y00",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        ); 

         volVectorField U
            (
               IOobject
            (
             "U",
             runTime.timeName(),
             mesh,
             IOobject::MUST_READ
            ), 
            mesh
            );



     // Done Reading Fields

     // Calculating Integral Quantities


        AInterface = 0;
        VolPhase0 = 0;
        VolPhase1 = 0;
        m00 = 0;
        m01 = 0;
        m10 = 0;
        m11 = 0;   
        velWaterPhase = zeroVec; //{0.0,0.0,0.0};
        velOilPhase = zeroVec; //{0.0,0.0,0.0};

        forAll(C10.internalField(), cellI)
        { 
            AInterface += A_intfc.internalField()[cellI];
            VolPhase0 += V[cellI]*(1 - alpha1.internalField()[cellI]);
            VolPhase1 += V[cellI]*alpha1.internalField()[cellI];
            m00 += V[cellI]*(1 - alpha1.internalField()[cellI])*rho0*Y00.internalField()[cellI];
            m01 += V[cellI]*(1 - alpha1.internalField()[cellI])*rho0*(1 - Y00.internalField()[cellI]);
            m10 += V[cellI]*alpha1.internalField()[cellI]*rho1*Y10.internalField()[cellI];
            m11 += V[cellI]*alpha1.internalField()[cellI]*rho1*(1 - Y10.internalField()[cellI]);
            velWaterPhase += V[cellI]*(1 - alpha1.internalField()[cellI])*U[cellI];
            velOilPhase += V[cellI]*alpha1.internalField()[cellI]*U[cellI];
        }                


           // Calculating outflux from reactor

        mDot00 = 0;
        mDot01 = 0;
        mDot10 = 0;
        mDot11 = 0;

       label frontBackID = mesh.boundaryMesh().findPatchID("frontAndBack");
         

        forAll(mesh.boundaryMesh(), outletID)
        {
            Info << "Patch " << outletID << endl;

            if(outletID != frontBackID)
            {

        forAll(mesh.boundaryMesh()[outletID], faceI)
        {
            mDot00 += (1 - alpha1.boundaryField()[outletID][faceI])*rho0*Y00.boundaryField()[outletID][faceI]*U.boundaryField()[outletID][faceI]&mesh.Sf().boundaryField()[outletID][faceI];
            mDot01 += (1 - alpha1.boundaryField()[outletID][faceI])*rho0*(1-Y00.boundaryField()[outletID][faceI])*U.boundaryField()[outletID][faceI]&mesh.Sf().boundaryField()[outletID][faceI];
            mDot10 += alpha1.boundaryField()[outletID][faceI]*rho1*Y10.boundaryField()[outletID][faceI]*U.boundaryField()[outletID][faceI]&mesh.Sf().boundaryField()[outletID][faceI];
            mDot11 += alpha1.boundaryField()[outletID][faceI]*rho1*(1 - Y10.boundaryField()[outletID][faceI])*U.boundaryField()[outletID][faceI]&mesh.Sf().boundaryField()[outletID][faceI];
        }

            }

            Info << "Done patch " << outletID << endl;

        }

           // Done Calculating outflux from reactor

     // Done Calculating Integral Quantities

     //  Calculating Vorticity Field


         volVectorField vorticity
            (
            IOobject
            (
                "vorticity",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            fvc::curl(U)
            );

        volScalarField magVorticity
            (
            IOobject
            (
                "magVorticity",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mag(vorticity)
            );

        magVorticity.write();

     // Done Calculating Vorticity Field


     // Writing to file ist)vs_t

        os<< runTime.value() << "    " << AInterface << "    " <<VolPhase0 << "    " << VolPhase1 << "    " << m00 << "    " << m01 << "    " << m10 << "    " << m11 << "    " << velWaterPhase.y() << "    " << velOilPhase.y() << "    " << mDot00 << "    " << mDot01  << "    " << mDot10  << "    " << mDot11  << "    "  << endl;

     // Done Writing to file ist)vs_t

        runTime.write();
    }        
   
    Info<< "ExecutionTime = "
        << runTime.elapsedCpuTime()
        << " s\n\n" << endl;     

    Info<< "End\n" << endl;    

    return 0;
}


// ************************************************************************* //
