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
#include "OFstream.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"            

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fileName outputFile("dropC10_vs_t");
    OFstream os(runTime.path()/outputFile);

    //scalar pi = Foam::constant::mathematical::pi;

    //scalar dz;
    //dz = runTime.controlDict().lookupOrDefault("dz",0.0001);    

    scalar ATot;
    scalar C10Tot;
    scalar C10Tot_0;
    scalar ATot_0;     

    const scalarField& V = mesh.V();    

    while (runTime.run())
    {        
        runTime++;        

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //Info<< "Reading alpha.oil field\n" << endl;        

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

        if((runTime.value() - 0.005) < 1e-5)
        {
            ATot_0 = 0;
            C10Tot_0 = 0;    
            forAll(C10.internalField(), cellI)
            {                 
                ATot_0 += A_intfc.internalField()[cellI];
                C10Tot_0 += V[cellI]*C10.internalField()[cellI];
            }
        }
               
        ATot = 0;
        C10Tot = 0;        
        
        forAll(C10.internalField(), cellI)
        { 
            C10Tot += V[cellI]*C10.internalField()[cellI];
            ATot += A_intfc.internalField()[cellI];                        
        }                
        
        os<< runTime.value() << "    " << C10Tot/C10Tot_0 << "    " << ATot/ATot_0
          << endl;

        runTime.write();
    }        
   
    Info<< "ExecutionTime = "
        << runTime.elapsedCpuTime()
        << " s\n\n" << endl;     

    Info<< "End\n" << endl;    

    return 0;
}


// ************************************************************************* //
