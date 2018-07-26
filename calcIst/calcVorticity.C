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

    fileName outputFile("ist_vs_t");
    OFstream os(runTime.path()/outputFile);

    //scalar pi = Foam::constant::mathematical::pi;

    //scalar dz;
    //dz = runTime.controlDict().lookupOrDefault("dz",0.0001);    

    scalar V0;
    scalar V1;
    scalar AInt;
    scalar m01;
    scalar m00;
    scalar m10;
    scalar m11;     

    const scalarField& V = mesh.V();    

    while (runTime.run())
    {        
        runTime++;        

        Info<< "Time = " << runTime.timeName() << nl << endl;


     IOobject Uheader
         (
             "U",
             runTime.timeName(),
             mesh,
             IOobject::MUST_READ
         );

     if (Uheader.headerOk())
     {
         Info<< "    Reading U" << endl;
         volVectorField U(Uheader, mesh);

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
        Info<< "Vorticity calculated and written"  << nl << endl;

     }
     else
     {
         Info<< "    No U" << endl;
     }
   
   

        runTime.write();
    }        
   
    Info<< "ExecutionTime = "
        << runTime.elapsedCpuTime()
        << " s\n\n" << endl;     

    Info<< "End\n" << endl;    

    return 0;
}


// ************************************************************************* //
