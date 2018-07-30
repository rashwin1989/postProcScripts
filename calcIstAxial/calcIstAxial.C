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
    calcIstAxial

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

    string outputFileName = "ist_Axial";
    string outputFileExt  = runTime.controlDict().lookup("caseName");
    outputFileName        = outputFileName + outputFileExt;
    fileName outputFile1(outputFileName);
    OFstream os(runTime.path()/outputFile1);

    const vector zeroVec(0,0,0);
    const double sqrt2            = 1.41421356;
    const scalar rho0             = 275;
    const scalar rho1             = 900;
    const scalarField& V          = mesh.V();    
    const vectorField& coordinate = mesh.C();    
 
    scalar dx;           // Currently not used
    scalar dy;
    scalar dz;           // Currently not used
    scalar yTop;
    scalar yBottom;
    scalar xLeft;
    scalar xRight;
    scalar tSteady;
    dx       = runTime.controlDict().lookupOrDefault("dx",0.00025);    
    dy       = runTime.controlDict().lookupOrDefault("dy",0.00025);    
    dz       = runTime.controlDict().lookupOrDefault("dz",0.001);    
    yTop     = runTime.controlDict().lookupOrDefault("yTop",0.240);    
    yBottom  = runTime.controlDict().lookupOrDefault("yBottom",0);    
    xLeft    = runTime.controlDict().lookupOrDefault("xLeft",-0.020);    
    xRight   = runTime.controlDict().lookupOrDefault("xRight",0.020);    
    tSteady  = runTime.controlDict().lookupOrDefault("tSteady",1);    
    
    scalar nCellsInY = (yTop - yBottom)/dy;
    if(nCellsInY - floor(nCellsInY) != 0)  // Want integer # of cells.Can perhaps be done in a better way!
    { 
           Info << "nCellsInY = " << nCellsInY << ". Resolution not entered consistently in controlDict." << endl;
    }
    label n     = floor(nCellsInY);
    label nTime = 0;

    double AInterface[n];
    double VolPhase0[n];
    double VolPhase1[n];
    double m00[n];
    double m01[n];
    double m10[n];
    double m11[n];
    double velYWaterPhase[n];
    double velYOilPhase[n];
    double MS00[n];
    for (int i=0; i < n; i++)  // Initializing dynamic Arrays
    {
        AInterface[i]     = 0;
        VolPhase0[i]      = 0;
        VolPhase1[i]      = 0;
        m00[i]            = 0;
        m01[i]            = 0;
        m10[i]            = 0;
        m11[i]            = 0;
        velYWaterPhase[i] = 0;
        velYOilPhase[i]   = 0;
        MS00[i]           = 0;
    }

    // Starting time loop
    while (runTime.run())
    {        
        runTime++;        
        scalar t=runTime.value(); 
        if (t >= tSteady)   // Average after steady state is reached
        {

            nTime = nTime +1;

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

         volScalarField mS00
             (
                IOobject
             (
             "mS00",
             runTime.timeName(),
             mesh,
             IOobject::MUST_READ
             ),
             mesh
             );

        // Done Reading Fields
        
        // Calculating Derived Fields

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
                IOobject::NO_WRITE
            ),
            mag(vorticity)
            );

         volScalarField QCriterion
            (
            IOobject
            (
                "QCriterion",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            0.5*mag(skew(fvc::grad(U)))*mag(skew(fvc::grad(U))) - mag(symm(fvc::grad(U)))*mag(symm(fvc::grad(U)))
            );

         volScalarField magStrainRate
            (
            IOobject
            (
                "magStrainRate",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sqrt2*mag(symm(fvc::grad(U)))
            );

            // Done Calculating Derived Fields

            // Calculating x-Averaged Quantities

            forAll(C10.internalField(), cellI)
            { 

                scalar xCoordinate = coordinate[cellI].x();
                scalar yCoordinate = coordinate[cellI].y(); 
                  /* Average only in main reactor, not inlet runners, etc. */
                if((xCoordinate>xLeft)&&(xCoordinate<xRight)&&(yCoordinate>yBottom)&&(yCoordinate<yTop)) 
                {
                    scalar yIndex = (coordinate[cellI].y() - 0)/dy;
                    /*   if(yIndex - floor(nCellsInY) != 0)
                    {
                        Info << "yIndex = " << yIndex << ". Resolution not entered consistently in controlDict." << endl;
                    } */
                    label yi = floor(yIndex);
                    /* Add quantities for all time steps, normalized after time loop ends */
                    AInterface[yi] = AInterface[yi] + A_intfc.internalField()[cellI];
                    VolPhase1[yi]  = VolPhase1[yi] + V[cellI]*alpha1.internalField()[cellI];
                    VolPhase0[yi]  = VolPhase0[yi] +  V[cellI]*(1 - alpha1.internalField()[cellI]);
                    m00[yi]  = m00[yi] + V[cellI]*(1 - alpha1.internalField()[cellI])*rho0*Y00.internalField()[cellI];
                    m01[yi]  = m01[yi] + V[cellI]*(1 - alpha1.internalField()[cellI])*rho0*(1 - Y00.internalField()[cellI]);
                    m10[yi]  = m10[yi] + V[cellI]*alpha1.internalField()[cellI]*rho1*Y10.internalField()[cellI];
                    m11[yi]  = m11[yi] + V[cellI]*alpha1.internalField()[cellI]*rho1*(1 - Y10.internalField()[cellI]);
                    /* Phase Velocity is weighted by volume of phase in that cross-section */
                    velYWaterPhase[yi]  =  velYWaterPhase[yi] + V[cellI]*(1 - alpha1.internalField()[cellI])*U[cellI].y();
                    velYOilPhase[yi]  = velYOilPhase[yi] + V[cellI]*alpha1.internalField()[cellI]*U[cellI].y(); 
                    MS00[yi]  = MS00[yi] + V[cellI]*mS00.internalField()[cellI];
                }
            }         // Done Calculating x-Averaged Quantities

        Info << "Time = " << runTime.timeName() << nl << endl;
        }   
    }   // Time loop over


    // Normalize quantities by # of time steps used for averaging
    for (int i=0; i < n; i++)
    {
        AInterface[i] = AInterface[i]/nTime;
        VolPhase0[i] = VolPhase0[i]/nTime;
        VolPhase1[i] = VolPhase1[i]/nTime;
        m00[i] = m00[i]/nTime;
        m01[i] = m01[i]/nTime;
        m10[i] = m10[i]/nTime;
        m11[i] = m11[i]/nTime;
        /* Phase Velocity is weighted by volume of phase in that cross-section */
        velYWaterPhase[i] = (velYWaterPhase[i]/nTime)/VolPhase0[i];
        velYOilPhase[i] = (velYOilPhase[i]/nTime)/VolPhase1[i];
        MS00[i] = MS00[i]/nTime;

        scalar yCoordinate = i*dy;
        os<< yCoordinate << "    " <<  AInterface[i] << "    " <<VolPhase0[i] << "    " << VolPhase1[i] << "    " << m00[i] << "    " << m01[i] << "    " << m10[i] << "    " << m11[i] << "    " << velYWaterPhase[i] << "    "  << velYOilPhase[i] << "    " << MS00[i]  << endl;
    }

   
    Info<< "ExecutionTime = "
        << runTime.elapsedCpuTime()
        << " s\n\n" << endl;     

    Info<< "End\n" << endl;    

    return 0;
}


// ************************************************************************* //
