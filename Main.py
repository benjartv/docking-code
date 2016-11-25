

from objects_py.TimeManager import *
from objects_py.ScoringFunction import *
from logic_py import EvolutionMethods
from utils_py import IO, Methods
import copy
import os
import random


'************************** PDB FILES ***************************'


__MoleculeName = '1ENY'
__LigandName = 'NAD'
#__testPDB = '1ENY_Rigid_RMSD_0.000.pdb'
#__testPDB = '1ENY_Rigid_RMSD_1.400.pdb'
__testPDB = '1ENY_Flexible_2_RMSD_1.300.pdb'
#__testPDB = '1ENY_Flexible_4_RMSD_1.400.pdb'
#__searchCenterPoint = [-4.200, 31.500, 12.000]
__searchCenterPoint = [-2.721, 33.882, 14.162] #ligand Center
__rotateAtoms = [[39,38], [35,36], [29,28], [28,27], [27,24], [24,23], [23,1], [1,4], [4,5], [5,6], [12,13]]
#[[4, 5], [5, 6]]#[[36, 35], [5, 6]] #[[23, 24], [36, 35], [12, 13], [5, 6]]

'''
__MoleculeName = '1AJV'
__LigandName = 'NMB'
#__testPDB = '1AJV_Rigid_RMSD_0.000.pdb'
#__testPDB = '1AJV_Rigid_RMSD_1.200.pdb'
#__testPDB = '1AJV_Flexible_2_RMSD_1.300.pdb'
__testPDB = '1AJV_Flexible_4_RMSD_1.400.pdb'
#__searchCenterPoint = [15.550, 21.000, 4.600]
__searchCenterPoint = [12.567, 22.694, 5.309] #ligand Center
__rotateAtoms = [[1,8], [3,23], [4,41], [5,30], [6,31], [7,32], [8,9], [16,17], [16,41], [23,24], [32,33], [33,34]]
#[[34, 33], [3, 23]] #[[34, 33], [16, 41], [8, 9], [3, 23]]


__MoleculeName = '1HPX'
__LigandName = 'KNI'
#__testPDB = '1HPX_Rigid_RMSD_0.000.pdb'
#__testPDB = '1HPX_Rigid_RMSD_1.100.pdb'
#__testPDB = '1HPX_Flexible_2_RMSD_1.200.pdb'
__testPDB = '1HPX_Flexible_4_RMSD_1.200.pdb'
#__searchCenterPoint = [4.000,  2.000, 10.500]
__searchCenterPoint = [6.099, 0.042, 12.521] #Ligand Center
__rotateAtoms = [[6,13],[10,11],[10,13],[15,16],[16,17],[16,20],[17,18],[22,25],[23,29],[25,30],[25,29],[26,27],[27,29],[37,40],[42,43]]
#[[26,27], [37,40]] #[[16,17], [26,27], [37,40], [42,43]]


__MoleculeName = '2UPJ'
__LigandName = 'U02'
#__testPDB = '2UPJ_Rigid_RMSD_0.000.pdb'
#__testPDB = '2UPJ_Rigid_RMSD_1.500.pdb'
#__testPDB = '2UPJ_Flexible_2_RMSD_0.900.pdb'
__testPDB = '2UPJ_Flexible_4_RMSD_1.300.pdb'
#__searchCenterPoint = [13.505, 23.500, 7.510]
__searchCenterPoint = [11.110, 24.008, 5.007] #Ligand Center
__rotateAtoms = [[4,5],[5,6],[8,9],[9,10],[10,11],[13,14],[18,20],[20,21],[20,24],[28,32],[30,31],[32,33],[32,35],[35,36]]
#[[14, 13], [36, 35]] #[[14,13], [20,24], [28,32], [36,35]]


__MoleculeName = '1BV9'
__LigandName = 'XV6'
#__testPDB = '1BV9_Rigid_RMSD_0.000.pdb'
#__testPDB = '1BV9_Rigid_RMSD_1.000.pdb'
#__testPDB = '1BV9_Flexible_2_RMSD_1.400.pdb'
__testPDB = '1BV9_Flexible_4_RMSD_1.700.pdb'
#__searchCenterPoint = [-8.000, 16.600, 24.600]
__searchCenterPoint = [-9.164, 15.912, 27.832] #Ligand Center
__rotateAtoms = [[2,9],[3,10],[4,7],[5,12],[7,27],[8,13],[11,20],[12,41],[13,14],[20,21],[29,33],[35,37],[43,47]]
#[[7,27],[41,12]] #[[7,27], [13,14], [41,12], [47,43]]
'''

'************************** Parameters **************************'

'Stopping criteria'
__generations = 300
__hours = 0

'Population'
__Population = 5
__searchSpaceSize = 5

'Algorithm'
__algorithmType = 2   # 1:Genetic Algorithm 2:Memetic Algorithm
__treeNodes = 3
__treeLevels = 3
__castPercent = [20, 50, 30]

'Mutation for Memetic'
__mutProbability = 0.1

'Distance criteria (for updateAgent)'
__distanceCriteria = 0.0

'Type of CrossOver'
__typeCO = 0
#0: Uniform CrossOver
#1: SPC CrossOver (window)
#2: 50/50 CrossOver (split and mix each block)
#3: Block CrossOver

'Memetic reset criteria'
__loopReset = 50

'Local Search'
__localSearch = True
__initTemperature = 1000
__minTemperature = 1
__alphaTemperature = 0.9
__innerLoop = 1
__typeLS = 0
#0: Normal LS
#1: LS with reduction by iteration
#2: LS with reduction by generation
#3: Mix of (2) and (3)

'Random Seed'
#random.seed(2139738276)


'****************************************************************'

'Timer start'
timer = TimeManager()
timer.start()

'Reads data'
moleculeList = IO.readAllPDB(__testPDB, __LigandName)
protein = copy.deepcopy(moleculeList[0])
ligand = copy.deepcopy(moleculeList[1])

'Center Molecule to (0,0,0)'
originPoint = protein.findCenter()
protein.translate([originPoint[0] * -1, originPoint[1] * -1, originPoint[2] * -1])
ligand.translate([originPoint[0] * -1, originPoint[1] * -1, originPoint[2] * -1])

'Traslate search center point to position of molecule'
alterCenterPoint = [__searchCenterPoint[0] + (originPoint[0]*-1),
                    __searchCenterPoint[1] + (originPoint[1]*-1),
                    __searchCenterPoint[2] + (originPoint[2]*-1)]

'Writes a new clean structure'
IO.writePDB("ligand.pdb", ligand, os.path.join(ROOTPATH,'temp'))
IO.writePDB("protein.pdb", protein, os.path.join(ROOTPATH,'temp'))

'Gets Connection Matrix'
connectMatrix = ligand.connectMatrix

'Init Scoring Function'
scoringFx = ScoringFunction()
scoringFx.defineLigandStructure("ligand.pdb", __LigandName, os.path.join(ROOTPATH, 'temp'))
ligand = IO.readPDB("ligand.pdb", 'HETATM', os.path.join(ROOTPATH, 'temp'))

'Recover Connection Matrix'
ligand.connectMatrix = connectMatrix


'Original'
IO.createsOriginalMolecule(__MoleculeName, __LigandName, connectMatrix)

'Generate Params'
__params = Methods.generatesParams(__searchSpaceSize,
                                   __searchCenterPoint,
                                   alterCenterPoint,
                                   __MoleculeName,
                                   __testPDB,
                                   __algorithmType,
                                   __treeNodes,
                                   __treeLevels,
                                   __castPercent,
                                   __rotateAtoms,
                                   __Population,
                                   __generations,
                                   __localSearch,
                                   __initTemperature,
                                   __minTemperature,
                                   __alphaTemperature,
                                   __innerLoop,
                                   __hours,
                                   __loopReset,
                                   __mutProbability,
                                   __typeLS,
                                   __typeCO,
                                   __distanceCriteria)


'Init Evolutionary algorithm'
EvolutionMethods.startAlgorithm(__params, scoringFx, protein, ligand)

'Timer Stop'
timer.stop()
timer.printTimeMinuts('*Timer: ')
