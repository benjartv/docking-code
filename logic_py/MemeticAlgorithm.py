

from objects_py.Gene import *
from objects_py.TimeManager import *
from logic_py import LocalSearch, ROOTPATH
from utils_py import Methods, IO
import rosetta
import random
import math
import copy
import os


# ------------------------ Memetic Algorithm ------------------------
class MemeticAlgorithm(object):

    # ------------------------------ Init -------------------------------
    def __init__(self, params, scoreFx, protein, ligand):
        self.__searchSpaceSize = params.searchSpaceSize
        self.__spaceCenter = params.searchCenterPoint
        self.__newSpaceCenter = params.newCenterPoint
        self.__rotateAtoms = params.rotatesAtoms
        self.__fileDir = params.folderName
        self.__moleculeName = params.moleculeName
        self.__ScoreFx = scoreFx
        self.__Protein = protein
        self.__Ligand = ligand
        self.__LocalSearch = None
        self.__isLocalSearch = params.isLocalSearch
        self.__generations = params.generationNumber
        self.__pocketSize = params.populationSize
        self.__loop = 0
        self.__log = ""
        self.__dataLog = ""
        self.__bestLoop = 0
        self.__loopReset = params.loopReset
        self.__lastReset = 0
        self.__resetCount = 0
        self.__treeNodes = params.populationCast[0]
        self.__treeLevels = params.populationCast[1]
        self.__treeSize = self.calculatesTree()
        self.__testPdb = params.testPdb
        self.__tmpDir = os.path.join(ROOTPATH,'temp')
        self.__treePopulation = [None] * (self.__treeSize * self.__pocketSize)
        self.__treeCurrent = [None] * self.__treeSize
        self.__timeStop = params.timeStop
        self.__angleList = params.angleList
        self.__mutProbability = params.mutProbability
        self.__typeLS = params.typeLS
        self.__typeCO = params.typeCO
        self.__distanceCriteria = params.distanceCriteria

        'Timer start'
        self.__timer = TimeManager()

        'Add Rotate Atoms'
        self.__Ligand.addRotationsAtoms(self.__rotateAtoms)

        'Pose Rosetta'
        self.__pose_ligand = self.__ScoreFx.getPoseFromPdbStr(self.__Protein, self.__Ligand)

        'Init Local Search'
        if(self.__isLocalSearch):
            self.__LocalSearch = LocalSearch.LocalSeach(params.initTemperature,
                                                        params.minTemperature,
                                                        params.alphaTemperature,
                                                        params.innerGenerations,
                                                        copy.copy(self.__ScoreFx),
                                                        self.__Protein,
                                                        self.__searchSpaceSize,
                                                        self.__newSpaceCenter,
                                                        self.__pose_ligand,
                                                        self.__typeLS,
							params.angleList)

    #--------------------------- Init Process ---------------------------
    def initProcess(self):
        self.__timer.start()
        self.initPopulation()

        while(self.validateEvolution()):
            self.generation()

            if(self.__ScoreFx.findBest()):
                self.__bestLoop = self.__loop


    #------------------------- Validates Iteration ----------------------
    def validateEvolution(self):

        self.printPopulation((self.__loop == 0 and "Init Population" or "Generation " + str(self.__loop)), True)

        self.__timer.stop()
        if(self.__generations > 0 and self.__loop >= self.__generations):
            self.generateFinalLog()
            return False
        elif(self.__timeStop > 0 and self.__timeStop <= self.__timer.getTimeHours()):
            self.generateFinalLog()
            return False
        else:
            self.__loop += 1
            self.reloadPopulation()
            return True


    #---------------------- Recalculate Population ----------------------
    def initPopulation(self):
        for i in xrange(self.__treeSize):
            cell = self.randomCell()
            cell = self.calculates(cell)
            self.addToCurrent(i,cell)
        self.updateTree(True)


    #---------------------- Recalculate Population ----------------------
    def generation(self):
        rootSize = int(self.__treeSize/self.__treeNodes)
        for i in xrange(rootSize):
            for j in xrange(self.__treeNodes):
                auxNodes = (i * self.__treeNodes) + (j + 1)
                p1 = self.randomSelection(i)
                p2 = self.randomSelection(auxNodes)
                if self.__typeCO == 1:
                    np = self.crossoverSPC(p1, p2)
                elif self.__typeCO == 2:
                    np = self.crossover50(p1, p2)
                elif self.__typeCO == 3:
                    np = self.crossoverBlock(p1, p2)
                else:
                    np = self.crossover(p1, p2)
                np = self.mutation(np)
                np = self.calculates(np)
                self.addToCurrent(auxNodes,np)
        self.updateTree()


    #------------------------ Using Local Seach ------------------------
    def calculates(self, np):
        if(self.__isLocalSearch):
            alphagen = (self.__generations - self.__loop) / float(self.__generations)
            np = self.__LocalSearch.initLocalProcess(np, self.__Ligand, alphagen)
            auxLigand = copy.deepcopy(self.__Ligand)
            auxLigand.traslateToPoint([(self.__newSpaceCenter[0] + np.x),
                                       (self.__newSpaceCenter[1] + np.y),
                                       (self.__newSpaceCenter[2] + np.z)])
            auxLigand.rotateByVector(Methods.spherePoint(1, np.sph_theta, np.sph_phi), np.theta)
            for r in xrange(len(self.__rotateAtoms)):
                auxLigand.rotateSegment(r, np.rotateAtoms[r])
            self.__ScoreFx.traceScoring(auxLigand,np.score)
        else:
            auxLigand = copy.deepcopy(self.__Ligand)
            auxLigand.traslateToPoint([(self.__newSpaceCenter[0] + np.x),
                                       (self.__newSpaceCenter[1] + np.y),
                                       (self.__newSpaceCenter[2] + np.z)])
            auxLigand.rotateByVector(Methods.spherePoint(1, np.sph_theta, np.sph_phi), np.theta)
            for r in xrange(len(self.__rotateAtoms)):
                auxLigand.rotateSegment(r, np.rotateAtoms[r])

            'Score'
            auxPose = rosetta.core.pose.Pose()
            auxPose.assign(self.__pose_ligand)
            for atom in xrange(auxPose.residue(auxPose.total_residue()).natoms()):
                atomV = rosetta.numeric.xyzVector_Real()
                atomV.x = round(auxLigand.x[atom],3)
                atomV.y = round(auxLigand.y[atom],3)
                atomV.z = round(auxLigand.z[atom],3)
                auxPose.residue(auxPose.total_residue()).set_xyz(auxLigand.atom[atom], atomV)
            np.score = self.__ScoreFx.generateScoringByPose(auxPose, auxLigand)
        return np

    #----------------------Generate pose for RMSD calculation------------

    def generatePose(self, np):
        auxLigand = copy.deepcopy(self.__Ligand)
        auxLigand.traslateToPoint([(self.__newSpaceCenter[0] + np.x),
                                   (self.__newSpaceCenter[1] + np.y),
                                   (self.__newSpaceCenter[2] + np.z)])
        auxLigand.rotateByVector(Methods.spherePoint(1, np.sph_theta, np.sph_phi), np.theta)
        for r in xrange(len(self.__rotateAtoms)):
            auxLigand.rotateSegment(r, np.rotateAtoms[r])

        'Score'
        auxPose = rosetta.core.pose.Pose()
        auxPose.assign(self.__pose_ligand)
        for atom in xrange(auxPose.residue(auxPose.total_residue()).natoms()):
            atomV = rosetta.numeric.xyzVector_Real()
            atomV.x = round(auxLigand.x[atom],3)
            atomV.y = round(auxLigand.y[atom],3)
            atomV.z = round(auxLigand.z[atom],3)
            auxPose.residue(auxPose.total_residue()).set_xyz(auxLigand.atom[atom], atomV)
        return auxPose 

    def distanceSolution(self, sol1, sol2):
        ligand1 = self.generatePose(sol1)
        ligand2 = self.generatePose(sol2)
        rmsd = rosetta.all_atom_rmsd(ligand1, ligand2)
        return rmsd

    #---------------------- Recalculate Population ----------------------
    def updateTree(self, initGeneration = False ):
        rootSize = int(self.__treeSize/self.__treeNodes)

        'Move currents to pocket'
        for a in xrange(len(self.__treeCurrent)):
            if(self.__treeCurrent[a] != None):
                self.currentToPocket(a)

        '''
        'Move best solution to root'
        if(not initGeneration):
            for i in xrange(rootSize - 1, -1, -1):
                for j in xrange(self.__treeNodes):
                    auxNode = (i * self.__treeNodes) + (j+1)
                    self.replaceBestNodes(i, auxNode)
        '''

        'Move best solution to root'
        if(not initGeneration):
            for i in xrange(rootSize - 1, -1, -1):
                node_ids = [(i * self.__treeNodes) + ( j + 1) for j in xrange(self.__treeNodes)]
                self.replaceFromLeafs(i, node_ids)


        'Recalculates Pockets'
        for j in xrange(self.__treeSize):
            auxPop = self.getPocketsFromNode(j)
            auxPop = self.recalculatesPocket(auxPop)
            self.setPocketsInNode(j,auxPop)


    #----------------


    #--------------------------- Add to Pocket --------------------------
    def currentToPocket(self, nodeId):
        newCell = copy.deepcopy(self.__treeCurrent[nodeId])
        self.__treeCurrent[nodeId] = None
        self.addToPocket(nodeId,newCell)


    #-------------------------- Add to Pocket ---------------------------
    '''
    def addToPocket(self, nodeId, newCell):
        population = []
        populationAux = self.getPocketsFromNode(nodeId)

        'Discart all white spaces'
        for x in xrange(len(populationAux)):
            if(populationAux[x]!=None):
                population.append(populationAux[x])

        'Replace '
        if(len(population) == self.__pocketSize ):
            population.sort(key=lambda x: x.score, reverse=False)
            for l in xrange(len(population)-1, -1, -1):
                if(population[l].score > newCell.score):
                    population[l] = newCell
                    break
        else:
            population.append(newCell)

        'Order Pocket'
        population.sort(key=lambda x: x.score, reverse=False)

        'Set Pockets in population'
        self.setPocketsInNode(nodeId, population)

    '''
    def addToPocket(self, nodeId, newCell):
        population = []
        populationAux = self.getPocketsFromNode(nodeId)

        'Discart all white spaces'
        for x in xrange(len(populationAux)):
            if(populationAux[x]!=None):
                population.append(populationAux[x])

        #(1) agent is empty
        if len(population) == 0:
            population.append(newCell)
        else:
            population.sort(key=lambda x: x.score, reverse=True)
            dcont = 0
            dArray = []
            #calculate distance between new solution and the pocket
            for l in xrange(0, len(population)):
                distance = self.distanceSolution(population[l], newCell)
                dArray.append(distance)
                if distance >= self.__distanceCriteria:
                    dcont+=1
            #agent is not complete
            if len(population) < self.__pocketSize:
                #(2) distance of new solution is higher than distance criteria for all solutions in pocket
                if dcont == len(population):
                    population.append(newCell)
                #(3) New solution is similar to the others, but the score is better than the best solution in the pocket
                elif ((dcont < len(population)) and (newCell.score < population[-1].score)):
                    population[dArray.index(min(dArray))] = newCell
            #agent is complete
            else:
                #(4) new solution score is better than the best of the pocket
                if newCell.score < population[-1].score:
                    population[dArray.index(min(dArray))] = newCell
                elif (newCell.score >= population[-1].score) and (newCell.score < population[0].score) and (dcont == len(population)):
                    #(5) distance  criteria is ok, replace the worst solution
                    population[0] = newCell

        'Order Pocket'
        population.sort(key=lambda x: x.score, reverse=False)

        'Set Pockets in population'
        self.setPocketsInNode(nodeId, population)


    #----------------------- Recalculates Pockets -----------------------
    def recalculatesPocket(self, population):
        popAux = []

        for q in xrange(len(population)):
            if(population[q] != None):
                popAux.append(copy.deepcopy(population[q]))

        popAux.sort(key=lambda x: x.score, reverse=False)

        return popAux[:]


    #------------------------- Replace Best Node ------------------------
    def replaceBestNodes(self, father, son):
        bestSonCell = son*self.__pocketSize
        bestFatherCell = father*self.__pocketSize

        if(self.__treePopulation[bestSonCell].score < self.__treePopulation[bestFatherCell].score):
            aux = copy.deepcopy(self.__treePopulation[bestSonCell])
            self.__treePopulation[bestSonCell] = self.__treePopulation[bestFatherCell]
            if(father == 0):
                self.addToPocket(father,aux)
            else:
                self.__treePopulation[bestFatherCell] = aux


    #------------------------ Replace from Leafs ------------------------
    def replaceFromLeafs(self, root_id, node_ids):

        bestFatherCell = root_id * self.__pocketSize
        bestSonCell = -1

        for i in xrange(len(node_ids)):
            auxNodeId = node_ids[i]*self.__pocketSize
            if(i == 0 or self.__treePopulation[auxNodeId].score < self.__treePopulation[bestSonCell].score):
                bestSonCell = auxNodeId

        if(bestSonCell != -1 and self.__treePopulation[bestSonCell].score < self.__treePopulation[bestFatherCell].score):
            aux = copy.deepcopy(self.__treePopulation[bestSonCell])
            self.__treePopulation[bestSonCell] = self.__treePopulation[bestFatherCell]
            if(root_id == 0):
                self.addToPocket(root_id, aux)
            else:
                self.__treePopulation[bestFatherCell] = aux


    #--------------------------- Add to Pocket --------------------------
    def addToCurrent(self, nodeId, newCell):
        self.__treeCurrent[nodeId] = copy.deepcopy(newCell)


    #------------------------- Gets Node Pockets ------------------------
    def getPocketsFromNode(self, nodeId):
        return self.__treePopulation[(nodeId*self.__pocketSize) : ((nodeId*self.__pocketSize)+self.__pocketSize)]


    #------------------------- Set Node Pockets -------------------------
    def setPocketsInNode(self, nodeId, pocket):
        for y in xrange(len(pocket)):
            self.__treePopulation[(nodeId*self.__pocketSize)+y] = copy.deepcopy(pocket[y])


    #-------------------------- Calculates Tree -------------------------
    def calculatesTree(self):
        totalNodes = 0
        for i in xrange(self.__treeLevels):
            totalNodes += self.__treeNodes**i
        return totalNodes


    # ------------------------- Apply Changes ---------------------------
    def applyChanges(self, solution):
        auxLigand = copy.deepcopy(self.__Ligand)
        auxLigand.traslateToPoint([(self.__newSpaceCenter[0] + solution.x),
                                   (self.__newSpaceCenter[1] + solution.y),
                                   (self.__newSpaceCenter[2] + solution.z)])
        auxLigand.rotateByVector(Methods.spherePoint(1, solution.sph_theta, solution.sph_phi), solution.theta)
        for r in xrange(len(self.__Ligand.rotationsPoints)):
            auxLigand.rotateSegment(r, solution.rotateAtoms[r])

        return auxLigand


    #------------------------- Random Selection --------------------------
    def randomSelection(self, nodeId):
        population = self.getPocketsFromNode(nodeId)
        popAux = []

        for q in xrange(len(population)):
            if(population[q] != None):
                popAux.append(copy.deepcopy(population[q]))

        randomId = random.randint(0, (len(popAux) - 1))
        selectedPop = copy.deepcopy(popAux[randomId])

        return selectedPop


    #---------------------- Recalculate Population ----------------------
    def randomCell(self):
        gen = Gene()
        gen.x = random.uniform(-self.__searchSpaceSize, self.__searchSpaceSize)
        gen.y = random.uniform(-self.__searchSpaceSize, self.__searchSpaceSize)
        gen.z = random.uniform(-self.__searchSpaceSize, self.__searchSpaceSize)
        gen.sph_theta = math.pi * random.randint(0, 200) / 100.0
        gen.sph_phi = math.pi * random.randint(0, 100) / 100.0
        gen.theta = math.pi * random.randint(0, 200) / 100.0

        for r in xrange(len(self.__rotateAtoms)):
            gen.rotateAtoms.append(math.pi * random.randint(0, 200) / 100.0)

        return copy.deepcopy(gen)


    #---------------------------- CrossOver -----------------------------
    def crossover(self, selectedPop1, selectedPop2 ):

        'Uniform Random Position'
        if(random.randint(0,1)==1):
            aux_x = selectedPop1.x
            selectedPop1.x = selectedPop2.x
            selectedPop2.x = aux_x

        if(random.randint(0,1)==1):
            aux_y = selectedPop1.y
            selectedPop1.y = selectedPop2.y
            selectedPop2.y = aux_y

        if(random.randint(0,1)==1):
            aux_z = selectedPop1.z
            selectedPop1.z = selectedPop2.z
            selectedPop2.z = aux_z

        'Uniform Random Vector'
        if(random.randint(0,1)==1):
            auxSphPhi = selectedPop1.sph_phi
            selectedPop1.sph_phi = selectedPop2.sph_phi
            selectedPop2.sph_phi = auxSphPhi

        if(random.randint(0,1)==1):
            auxSphTheta = selectedPop1.sph_theta
            selectedPop1.sph_theta = selectedPop2.sph_theta
            selectedPop2.sph_theta= auxSphTheta

        'Switch Theta'
        if(random.randint(0,1)==1):
            aux_theta = selectedPop1.theta
            selectedPop1.theta = selectedPop2.theta
            selectedPop2.theta = aux_theta

        'Switch Bonds'
        rot_bond0 = len(selectedPop1.rotateAtoms)
        rot_bond1 = len(selectedPop2.rotateAtoms)
        if( rot_bond0 > 0 and rot_bond1 > 0 and rot_bond0 == rot_bond1):
            for i in xrange(rot_bond0):
                if(random.randint(0,1)==1):
                    aux_bond = selectedPop1.rotateAtoms[i]
                    selectedPop1.rotateAtoms[i] = selectedPop2.rotateAtoms[i]
                    selectedPop2.rotateAtoms[i] = aux_bond

        return selectedPop1

    #---------------------------- CrossOver Block-----------------------------
    # split the solution in 3 parts: center of the ligand, rotational parameters and 
    # dihedral angles (all of them)
    def crossoverBlock(self, selectedPop1, selectedPop2):
        newPop = Gene()
        'Switch center of the ligand'
        if(random.randint(0,1) == 1):
            newPop.x = selectedPop1.x
            newPop.y = selectedPop1.y
            newPop.z = selectedPop1.z
        else:
            newPop.x = selectedPop2.x
            newPop.y = selectedPop2.y
            newPop.z = selectedPop2.z

        'Switch sphere coordinates'
        if(random.randint(0,1) == 1):

            newPop.sph_theta = selectedPop1.sph_theta
            newPop.sph_phi = selectedPop1.sph_phi
            newPop.theta = selectedPop1.theta
        else:
            newPop.sph_theta = selectedPop2.sph_theta
            newPop.sph_phi = selectedPop2.sph_phi
            newPop.theta = selectedPop2.theta

        'Switch dihedral bonds'
        if(random.randint(0,1) == 1):
            newPop.rotateAtoms = selectedPop1.rotateAtoms[:]
        else:
            newPop.rotateAtoms = selectedPop2.rotateAtoms[:]

        return copy.deepcopy(newPop)
        '''
    def crossoverBlock(self, selectedPop1, selectedPop2):
        'Switch center of the ligand'
        if(random.randint(0,1) == 1):
            aux_x = selectedPop1.x
            aux_y = selectedPop1.y
            aux_z = selectedPop1.z
            selectedPop1.x = selectedPop2.x
            selectedPop1.y = selectedPop2.y
            selectedPop1.z = selectedPop2.z
            selectedPop2.x = aux_x
            selectedPop2.y = aux_y
            selectedPop2.z = aux_z

        'Switch sphere coordinates'
        if(random.randint(0,1) == 1):
            aux_sph_phi = selectedPop1.sph_phi
            aux_sph_theta = selectedPop1.sph_theta
            aux_theta = selectedPop1.theta
            selectedPop1.sph_phi = selectedPop2.sph_phi
            selectedPop1.sph_theta = selectedPop2.sph_theta
            selectedPop1.theta = selectedPop2.theta
            selectedPop2.sph_phi = aux_sph_phi
            selectedPop2.sph_theta = aux_sph_theta
            selectedPop2.theta = aux_theta

        'Switch dihedral bonds'
        if(random.randint(0,1) == 1):
            aux_bonds = selectedPop1.rotateAtoms[:]
            selectedPop1.rotateAtoms = selectedPop2.rotateAtoms[:]
            selectedPop2.rotateAtoms = aux_bonds[:]

        return selectedPop1
        '''
    #---------------------------- CrossOver SPC-----------------------------
    # window of gene and mix the solutions
    def crossoverSPC(self, selectedPop1, selectedPop2):

        popSolution1 = []
        popSolution2 = []

        popSolution1.append(selectedPop1.x)
        popSolution1.append(selectedPop1.y)
        popSolution1.append(selectedPop1.z)
        popSolution1.append(selectedPop1.sph_theta)
        popSolution1.append(selectedPop1.sph_phi)
        popSolution1.append(selectedPop1.theta)

        for angle in selectedPop1.rotateAtoms:
            popSolution1.append(angle)

        popSolution2.append(selectedPop2.x)
        popSolution2.append(selectedPop2.y)
        popSolution2.append(selectedPop2.z)
        popSolution2.append(selectedPop2.sph_theta)
        popSolution2.append(selectedPop2.sph_phi)
        popSolution2.append(selectedPop2.theta)

        for angle in selectedPop2.rotateAtoms:
            popSolution2.append(angle)

        genSize = len(popSolution1)
        cutPoint = random.randint(0, genSize-1)
        length = random.randint(0, genSize-1)
        if length == 0:
            return selectedPop1
        elif length > (genSize - cutPoint) :
            newpopSolution = popSolution1[cutPoint:]
            turnTop = length-(genSize-cutPoint)
            for i in range(0,turnTop):
                newpopSolution.insert(i,popSolution1[i])
            for i in range(turnTop, cutPoint):
                newpopSolution.insert(i,popSolution2[i])
        else:
            newpopSolution = popSolution1[cutPoint:cutPoint+length]
            for i in range(0,cutPoint):
                newpopSolution.insert(i,popSolution2[i])
            for i in range(cutPoint+length,genSize):
                newpopSolution.insert(i,popSolution2[i])

        newPop = Gene()
        newPop.x = newpopSolution[0]
        newPop.y = newpopSolution[1]
        newPop.z = newpopSolution[2]
        newPop.sph_theta = newpopSolution[3]
        newPop.sph_phi = newpopSolution[4]
        newPop.theta = newpopSolution[5]
        newPop.rotateAtoms = newpopSolution[6:]
        return copy.deepcopy(newPop)

    '''
    def crossoverSPC(self, selectedPop1, selectedPop2):

        popSolution1 = []
        popSolution2 = []

        popSolution1.append(selectedPop1.x)
        popSolution1.append(selectedPop1.y)
        popSolution1.append(selectedPop1.z)
        popSolution1.append(selectedPop1.sph_phi)
        popSolution1.append(selectedPop1.sph_theta)
        popSolution1.append(selectedPop1.theta)

        for angle in selectedPop1.rotateAtoms:
            popSolution1.append(angle)

        popSolution2.append(selectedPop2.x)
        popSolution2.append(selectedPop2.y)
        popSolution2.append(selectedPop2.z)
        popSolution2.append(selectedPop2.sph_phi)
        popSolution2.append(selectedPop2.sph_theta)
        popSolution2.append(selectedPop2.theta)

        for angle in selectedPop2.rotateAtoms:
            popSolution2.append(angle)

        genSize = len(popSolution1)
        cutPoint = random.randint(0, genSize-1)
        length = random.randint(0, genSize-1)
        if length == 0:
            return selectedPop1
        elif length > (genSize - cutPoint) :
            newpopSolution = popSolution1[cutPoint:]
            turnTop = length-(genSize-cutPoint)
            for i in range(0,turnTop):
                newpopSolution.insert(i,popSolution1[i])
            for i in range(turnTop, cutPoint):
                newpopSolution.insert(i,popSolution2[i])
        else:
            newpopSolution = popSolution1[cutPoint:cutPoint+length]
            for i in range(0,cutPoint):
                newpopSolution.insert(i,popSolution2[i])
            for i in range(cutPoint+length,genSize):
                newpopSolution.insert(i,popSolution2[i])

        selectedPop1.x = newpopSolution[0]
        selectedPop1.y = newpopSolution[1]
        selectedPop1.z = newpopSolution[2]
        selectedPop1.sph_phi = newpopSolution[3]
        selectedPop1.sph_theta = newpopSolution[4]
        selectedPop1.theta = newpopSolution[5]
        newbonds = newpopSolution[6:]
        for i in range(0,len(newbonds)):
            selectedPop1.rotateAtoms[i] = newbonds[i]

        return selectedPop1
    '''
    #---------------------------- CrossOver 50/50-----------------------------
    # from each block (center, rot, bonds) split the blocks and switch in to the parents
    def crossover50(self, selectedPop1, selectedPop2):
        'Switch center of the ligand'
        newPop = Gene()
        centerrand = random.randint(0,5)
        if centerrand == 0:
            newPop.x = selectedPop1.x
            newPop.y = selectedPop2.y
            newPop.z = selectedPop2.z
        elif centerrand == 1:
            newPop.x = selectedPop1.x
            newPop.y = selectedPop1.y
            newPop.z = selectedPop2.z
        elif centerrand == 2:
            newPop.x = selectedPop1.x
            newPop.y = selectedPop1.y
            newPop.z = selectedPop1.z
        elif centerrand == 3:
            newPop.x = selectedPop2.x
            newPop.y = selectedPop2.y
            newPop.z = selectedPop2.z
        elif centerrand == 4:
            newPop.x = selectedPop2.x
            newPop.y = selectedPop1.y
            newPop.z = selectedPop1.z
        else:
            newPop.x = selectedPop2.x
            newPop.y = selectedPop2.y
            newPop.z = selectedPop1.z
        'Switch rotation of the ligand'
        rotrand = random.randint(0,5)
        if rotrand == 0:
            newPop.sph_theta = selectedPop1.sph_theta
            newPop.sph_phi = selectedPop2.sph_phi
            newPop.theta = selectedPop2.theta
        elif rotrand == 1:
            newPop.sph_theta = selectedPop1.sph_theta
            newPop.sph_phi = selectedPop1.sph_phi
            newPop.theta = selectedPop2.theta
        elif rotrand == 2:
            newPop.sph_theta = selectedPop1.sph_theta
            newPop.sph_phi = selectedPop1.sph_phi
            newPop.theta = selectedPop1.theta
        elif rotrand == 3:
            newPop.sph_theta = selectedPop2.sph_theta
            newPop.sph_phi = selectedPop2.sph_phi
            newPop.theta = selectedPop2.theta
        elif rotrand == 4:
            newPop.sph_theta = selectedPop2.sph_theta
            newPop.sph_phi = selectedPop1.sph_phi
            newPop.theta = selectedPop1.theta
        else:
            newPop.sph_theta = selectedPop2.sph_theta
            newPop.sph_phi = selectedPop2.sph_phi
            newPop.theta = selectedPop1.theta
        'Switch dihedral bonds of the ligand'
        bondrand = random.randint(0,len(selectedPop1.rotateAtoms)-1)
        newbonds1 = selectedPop1.rotateAtoms[:bondrand] + selectedPop2.rotateAtoms[bondrand:]
        newbonds2 = selectedPop2.rotateAtoms[:bondrand] + selectedPop1.rotateAtoms[bondrand:]
        if(random.randint(0,1) == 1):
            newPop.rotateAtoms = newbonds1
        else:
            newPop.rotateAtoms = newbonds2

        return copy.deepcopy(newPop)

        '''
    def crossover50(self, selectedPop1, selectedPop2):
        'Switch center of the ligand'
        centerrand = random.randint(0,4)
        if centerrand == 0:
            aux_x = selectedPop1.x
            selectedPop1.x = selectedPop2.x
            selectedPop2.x = aux_x
        elif centerrand == 1:
            aux_x = selectedPop1.x
            aux_y = selectedPop1.y
            selectedPop1.x = selectedPop2.x
            selectedPop1.y = selectedPop2.y
            selectedPop2.x = aux_x
            selectedPop2.y = aux_y
        elif centerrand == 2:
            aux_y = selectedPop1.y
            aux_z = selectedPop1.z
            selectedPop1.y = selectedPop2.y
            selectedPop1.z = selectedPop2.z
            selectedPop2.y = aux_y
            selectedPop2.z = aux_z
        elif centerrand == 3:
            aux_z = selectedPop1.z
            selectedPop1.z = selectedPop2.z
            selectedPop2.z = aux_z
        else:
            aux_x = selectedPop1.x
            aux_y = selectedPop1.y
            aux_z = selectedPop1.z
            selectedPop1.x = selectedPop2.x
            selectedPop1.y = selectedPop2.y
            selectedPop1.z = selectedPop2.z
            selectedPop2.x = aux_x
            selectedPop2.y = aux_y
            selectedPop2.z = aux_z
        'Switch rotation of the ligand'
        rotrand = random.randint(0,4)
        if rotrand == 0:
            aux_sph_theta = selectedPop1.sph_theta
            selectedPop1.sph_theta = selectedPop2.sph_theta
            selectedPop2.sph_theta = aux_sph_theta
        elif rotrand == 1:
            aux_sph_theta = selectedPop1.sph_theta
            aux_sph_phi = selectedPop1.sph_phi
            selectedPop1.sph_theta = selectedPop2.sph_theta
            selectedPop1.sph_phi = selectedPop2.sph_phi
            selectedPop2.sph_theta = aux_sph_theta
            selectedPop2.sph_phi = aux_sph_phi
        elif rotrand == 2:
            aux_sph_phi = selectedPop1.sph_phi
            aux_theta = selectedPop1.theta
            selectedPop1.sph_phi = selectedPop2.sph_phi
            selectedPop1.theta = selectedPop2.theta
            selectedPop2.sph_phi = aux_sph_phi
            selectedPop2.theta = aux_theta
        elif rotrand == 3:
            aux_theta = selectedPop1.theta
            selectedPop1.theta = selectedPop2.theta
            selectedPop2.theta = aux_theta
        else:
            aux_sph_theta = selectedPop1.sph_theta
            aux_sph_phi = selectedPop1.sph_phi
            aux_theta = selectedPop1.theta
            selectedPop1.sph_theta = selectedPop2.sph_theta
            selectedPop1.sph_phi = selectedPop2.sph_phi
            selectedPop1.theta = selectedPop2.theta
            selectedPop2.sph_theta = aux_sph_theta
            selectedPop2.sph_phi = aux_sph_phi
            selectedPop2.theta = aux_theta
        'Switch dihedral bonds of the ligand'
        bondrand = random.randint(0,len(selectedPop1.rotateAtoms)-1)
        newbonds1 = selectedPop1.rotateAtoms[:bondrand] + selectedPop2.rotateAtoms[bondrand:]
        newbonds2 = selectedPop2.rotateAtoms[:bondrand] + selectedPop1.rotateAtoms[bondrand:]
        if(random.randint(0,1) == 1):
            selectedPop1.rotateAtoms = newbonds1[:]
            selectedPop2.rotateAtoms = newbonds2[:]
        else:
            selectedPop2.rotateAtoms = newbonds1[:]
            selectedPop1.rotateAtoms = newbonds2[:]

        return selectedPop1
        '''
    #----------------------------- Mutation -----------------------------
    def mutation(self, selectedPop):

        muationProb = self.__mutProbability

        if(random.uniform(0,1) <= muationProb):
            rotation_bonds = len(selectedPop.rotateAtoms)
            if(rotation_bonds > 0):
                mutPos = random.randint(1, (6 + rotation_bonds))
            else:
                 mutPos = random.randint(1, 6)
            if mutPos == 1:
                selectedPop.x = random.uniform(-self.__searchSpaceSize, self.__searchSpaceSize)
            elif mutPos == 2:
                selectedPop.y = random.uniform(-self.__searchSpaceSize, self.__searchSpaceSize)
            elif mutPos == 3:
                selectedPop.z = random.uniform(-self.__searchSpaceSize, self.__searchSpaceSize)
            elif mutPos == 4:
                selectedPop.sph_theta = math.pi * random.randint(0, 200) / 100.0
            elif mutPos == 5:
                selectedPop.sph_phi = math.pi * random.randint(0, 100) / 100.0
            elif mutPos == 6:
                selectedPop.theta = math.pi * random.randint(0, 200) / 100.0
            elif mutPos > 6:
                selectedPop.rotateAtoms[mutPos-7] = math.pi * random.randint(0, 200) / 100.0

        return selectedPop


    #------------------------- Print Population -------------------------
    def printPopulation(self, title, onLog = False):
        textList = []
        dataList = []

        textList.append("\n* ")
        textList.append(title)
        textList.append("\n")
        textList.append("**************************************************************************")
        for nodeId in xrange(self.__treeSize):
            population = self.getPocketsFromNode(nodeId)
            textList.append("\n")
            textList.append("Nodo ")
            textList.append(str(nodeId))
            textList.append("\n")

            for i in xrange(len(population)):
                textList.append(str(i+1))
                textList.append("-\t")
                textList.append("[\t")

                if(population[i] == None):
                    textList.append("None")
                    textList.append(" ]\n")
                    dataList.append("0")
                else:
                    textList.append("{0:.3f}".format(self.__newSpaceCenter[0] + population[i].x))
                    textList.append("\t")
                    textList.append("{0:.3f}".format(self.__newSpaceCenter[1] + population[i].y))
                    textList.append("\t")
                    textList.append("{0:.3f}".format(self.__newSpaceCenter[2] + population[i].z))
                    textList.append("\t")
                    textList.append("{0:.3f}".format(population[i].sph_theta))
                    textList.append("\t")
                    textList.append("{0:.3f}".format(population[i].sph_phi))
                    textList.append("\t")
                    textList.append("{0:.3f}".format(population[i].theta))
                    textList.append("\t")

                    for j in xrange(len(self.__rotateAtoms)):
                        textList.append("{0:.3f}".format(population[i].rotateAtoms[j]))
                        textList.append("\t")

                    textList.append(" ]\t")
                    textList.append("score: ")
                    textList.append("{0:.3f}".format(population[i].score))
                    textList.append("\t \t")
                    textList.append("\n")
                    dataList.append("{0:.3f}".format(population[i].score))

                if(nodeId != (self.__treeSize - 1) or i < (len(population) - 1)):
                    dataList.append(",")

        dataList.append("\n")


        if (onLog == False):
            print "\n* " + title + "\n"
            print ''.join(textList)
        else:
            self.__log += ''.join(textList)
            self.__dataLog += ''.join(dataList)


    #------------------------- Print Final Log --------------------------
    def generateFinalLog(self):
        self.__timer.stop()

        'Writes Files'
        ligand_PDB_name = "Ligand_score_" + "{0:.3f}".format(self.__ScoreFx.getBestScore()) + ".pdb"
        ligandProtein_PDB_name = "LigandProtein_score_" + "{0:.3f}".format(self.__ScoreFx.getBestScore()) + ".pdb"
        IO.writePDB(ligand_PDB_name, self.__ScoreFx.getBestLigand(), self.__fileDir)
        IO.writeAllPDB(ligandProtein_PDB_name, self.__Protein, self.__ScoreFx.getBestLigand(), self.__fileDir)

        initLog = "Molecule Name: " + self.__moleculeName + "\n"
        initLog += "Test Pdb: " + self.__testPdb + "\n"
        initLog += "Pocket: " + str(self.__pocketSize) + "\n"
        initLog += "Generations: " + str(self.__generations) + "\n"
        initLog += "Time to Stop: " + str(self.__timeStop) + "\n"
        initLog += "Space Size: " + str(self.__searchSpaceSize) + "\n"
        initLog += "Point: ( " + str(self.__spaceCenter[0]) + " , " + str(self.__spaceCenter[1]) + " , " + str(self.__spaceCenter[2]) + " )\n"
        initLog += "Local Search: " + str(self.__isLocalSearch) + "\n"
        initLog += "Type of Crossover: " + str(self.__typeCO) + "\n"
        initLog += "Type of Local Search: " + str(self.__typeLS) + "\n"
        initLog += "Mutation Probability: " + str(self.__mutProbability) + "\n"
        initLog += "Distance Criteria: " + str(self.__distanceCriteria) + "\n"
        initLog += "Rotate Atoms: " + str(self.__rotateAtoms) + "\n"
        initLog += "\n- Best Score: " + str(self.__ScoreFx.getBestScore()) + "\n"
        initLog += "- Reach Best Score: " + str(self.__bestLoop) + "\n"
        auxRMSD = str(self.__ScoreFx.getRMSD("original_" + self.__moleculeName + ".pdb", ligandProtein_PDB_name, self.__tmpDir, self.__fileDir))
        initLog += "- RMSD: " + auxRMSD + "\n"
        initLog += "- Timer: " + str(self.__timer.getTimeChronometer()) + "\n" + "\n"
        self.__log = initLog + self.__log

        'Data Log'
        IO.writeLog("iterations.log", self.__log, self.__fileDir)
        IO.writeLog("data.log", self.__dataLog, self.__fileDir)
        IO.writeHistoricLog("historic.log", str(self.__ScoreFx.getBestScore()) + ";" + auxRMSD + ";" + str(self.__bestLoop) + ";" + str(self.__loop) , self.__fileDir)

        'Delete files'
        Methods.deleteTempFiles(self.__ScoreFx.getLigandName(), self.__moleculeName)


    #----------------------- Reload Population --------------------------
    def reloadPopulation(self):
        if(self.__loopReset != 0):
            if(self.__lastReset >= self.__loopReset):
                self.__lastReset = 0
                self.__loopReset = self.__resetCount
                self.__resetCount = 0
                for i in xrange(1,len(self.__treePopulation)):
                    if( i % self.__pocketSize != 0):
                        self.__treePopulation[i] = None
            elif(self.__bestLoop == self.__loop):
                self.__lastReset = 0
                self.__resetCount += 1
            else:
                self.__lastReset += 1
                self.__resetCount += 1
