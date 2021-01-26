'''
Topic: Applying imperfections in form of buckling mode or random imperfection
By: Matheus C. Fernandes
Find latest version of this code at: http://fer.me/git/abaqus-scripting

'''
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import numpy as np
import random
random.seed(1111)

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)


def ApplyBuckling(mdbName,odbName,ImpFrames,ImpWeights,StepName,PartName,ModelName,InstanceName,AllNodes):
    ########## DESCRIPTION OF VARIABLES FOR THE DEFINED FUNCTION
    '''
    mdbName -- 	THIS IS THE NAME FOR THE MDB PATH FOR THE CAE
    odbName -- THIS IS THE PATH FOR THE ODB FILE
    ImpFrames -- THIS IS THE FRAME NUMBERS OF THE FREQUENCY NUMBER OF THE ANALYSIS
    ImpWeigts -- THIS IS THE WEIGHT FOR THE FREQUENCY DISPLACEMENT B.C.
    StepName -- THIS IS THE NAME OF THE STEP
    PartName -- THIS IS THE NAME OF THE PART WHERE THE BUCKLING WILL BE APPLIED
    ModelName -- THIS IS THE MODEL NAME IN THE CAE
    InstanceName -- THIS IS THE INSTANCE NAME IN THE CAE
    AllNodes -- THIS IS THE SET NAME THAT CONTAINS ALL OF THE NODES
    '''
    ##########
    # BEGIN CODE BY IMPORTING THE MDB MODEL OF FREQUENCY ANALYSIS
    mdb=openMdb(pathName=mdbName)
    odb=openOdb(path=odbName)
    # THIS OPENS THE RESULTS ODB WITH THE NODESET OF ALL OF THE NODES
    pbpPartNodes=odb.rootAssembly.instances[InstanceName.upper()].nodeSets[AllNodes.upper()]

    #CREATE A MATRIX TO SAVE THE NEW COORDINATES OF ALL NODES
    NewCoord=np.zeros((len(mdb.models[ModelName].parts[PartName].nodes), 3))	
    
    for CImp in range(len(ImpFrames)):
        cframe = ImpFrames[CImp]
        firstFrame = odb.steps[StepName].frames[cframe]
        displacement = firstFrame.fieldOutputs['U']
        pbpDispField = displacement.getSubset(region=pbpPartNodes)
        pbpDisp = pbpDispField.values


        # Imperfection Using Buckling Analysis Results
        #---------------------------------------------------------------
        ind=0;
        IMP = ImpWeights[CImp]
        #CALCULATE THE MODIFIED COORDINATES
        for i in mdb.models[ModelName].parts[PartName].nodes:
            NewCoord[ind][0]=i.coordinates[0]+IMP*pbpDisp[ind].data[0]
            NewCoord[ind][1]=i.coordinates[1]+IMP*pbpDisp[ind].data[1]
            #NewCoord[ind][2]=i.coordinates[2]+IMP*pbpDisp[ind].data[2]
            ind=ind+1

        #SET THE NEW COORDINATES
        mdb.models[ModelName].parts[PartName].editNode(
            nodes=mdb.models[ModelName].parts[PartName].nodes,
            coordinates=NewCoord)

    mdb.models[ModelName].rootAssembly.regenerate()

    print 'Original Coordinates Modified Successfully!'
 

def ApplyRandomImperfection(mdb,IMP,PartName,ModelName):
    ########## DESCRIPTION OF VARIABLES FOR THE DEFINED FUNCTION
    '''
    IMP -- THIS IS THE WEIGHT FOR THE FREQUENCY DISPLACEMENT B.C.
    StepName -- THIS IS THE NAME OF THE STEP
    PartName -- THIS IS THE NAME OF THE PART WHERE THE BUCKLING WILL BE APPLIED
    ModelName -- THIS IS THE MODEL NAME IN THE CAE
    InstanceName -- THIS IS THE INSTANCE NAME IN THE CAE
    AllNodes -- THIS IS THE SET NAME THAT CONTAINS ALL OF THE NODES
    '''
    ##########
    # BEGIN CODE BY IMPORTING THE MDB MODEL OF FREQUENCY ANALYSIS

    #CREATE A MATRIX TO SAVE THE NEW COORDINATES OF ALL NODES
    NewCoord=np.zeros((len(mdb.models[ModelName].parts[PartName].nodes), 3))	
    
    # Imperfection Using Buckling Analysis Results
    #---------------------------------------------------------------
    ind=0
    #CALCULATE THE MODIFIED COORDINATES
    for i in mdb.models[ModelName].parts[PartName].nodes:
        NewCoord[ind][0]=i.coordinates[0]+IMP*(0.5-random.random())
        NewCoord[ind][1]=i.coordinates[1]+IMP*(0.5-random.random())
        NewCoord[ind][2]=i.coordinates[2]#+IMP*(0.5-random.random())
        ind=ind+1

    #SET THE NEW COORDINATES
    mdb.models[ModelName].parts[PartName].editNode(
        nodes=mdb.models[ModelName].parts[PartName].nodes,
        coordinates=NewCoord)

    mdb.models[ModelName].rootAssembly.regenerate()

    print 'Original Coordinates Modified Successfully!'
