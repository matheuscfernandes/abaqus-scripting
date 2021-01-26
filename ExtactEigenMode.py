'''
Topic: Extracting eigenmode from simulations
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


def ExtractEigenMode(JobName, NumberOfModes):
    odb = openOdb(path=JobName + '.odb')
    Freq = []
    for i in xrange(1, NumberOfModes + 1):
        Desc = odb.getFrame(i).description
        Desc = Desc.split("=")
        Freq.append(float(Desc[1]))
    odb.close()
    return Freq


