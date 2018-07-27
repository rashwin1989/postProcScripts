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
    extractScalarFields

Description
    Extracts scalar fields along a line and writes to file as:
    y co-ordinate    field value
    .
    .
    .
    Should sort both vectors in Matlab before plotting

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "fileName.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"            

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    wordList fieldNames(runTime.controlDict().lookup("fieldsToExtract"));
    label n = fieldNames.size();
    Info<< "Number of fields to extract:  " << n << endl;

    word fixCoName(runTime.controlDict().lookup("fixedCoordinateName"));
    scalarList fixCoValues(runTime.controlDict().lookup("fixedCoordinateValues"));
    scalar dxValue; scalar dyValue; scalar dzValue; scalar zSliceValue;
    dxValue = runTime.controlDict().lookupOrDefault("dxValue",0.0005);
    dyValue = runTime.controlDict().lookupOrDefault("dyValue",0.00005);
    dzValue = runTime.controlDict().lookupOrDefault("dzValue",0.00005);
    zSliceValue = runTime.controlDict().lookupOrDefault("zSliceValue",0);
    label nFixCo = fixCoValues.size();
    Info<< "Number of " << fixCoName << " values to extract at:  " << nFixCo << endl;
    
    while (runTime.run())
    {           
        runTime++;        

        Info<< "Time = " << runTime.timeName() << nl << endl;

        for(label i=0; i<n; i++)
        {
            word curFieldName = fieldNames[i];
            Info<< "Reading field  " << curFieldName << endl;

            volScalarField curField
            (
                IOobject
                (
                    curFieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );            

            const scalarField& curFieldCells = curField.internalField();
            const vectorField& C = mesh.C();

            for(label j=0; j<nFixCo; j++)
            {
                scalar curFixCoValue = fixCoValues[j];
                Info<< "Extracting field at " << fixCoName << "= " << curFixCoValue << endl;

                fileName curOutputFile(curFieldName+"_"+fixCoName+"="+Foam::name(curFixCoValue));
                OFstream os(runTime.path()/runTime.timeName()/curOutputFile);
            
                forAll(curFieldCells,cellI)
                {
                    scalar xCell = C[cellI].x();
                    scalar yCell = C[cellI].y();
                    scalar zCell = C[cellI].z();                
                    if(fixCoName == "x")
                    {                    
                        if(mag(xCell - curFixCoValue) < 0.5*dxValue && mag(zCell - zSliceValue) < 0.5*dzValue)
                        {
                            os<< yCell << "    " << curFieldCells[cellI] << endl;
                        }
                    }
                    else
                    {
                        if(mag(yCell - curFixCoValue) < 0.5*dyValue && mag(zCell - zSliceValue) < 0.5*dzValue)
                        {
                            os<< xCell << "    " << curFieldCells[cellI] << endl;
                        }
                    }
                }
            }
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
