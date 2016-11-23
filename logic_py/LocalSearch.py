

from utils_py import Methods
import rosetta
import random
import math
import copy


# -------------------------- Local Search ---------------------------
class LocalSeach(object):

    # ------------------------------ Init -------------------------------
    def __init__(self, temp, tempMin, tempAlpha, innerLoop,  scoreFx, protein, spaceLen, spaceCenter, poseLigand, typeLS, angleList):
        self.__temp = temp
        self.__tempMin = tempMin
        self.__tempAlpha = tempAlpha
        self.__innerIteration = innerLoop
        self.__scoreFx = scoreFx
        self.__protein = protein
        self.__spaceLen = spaceLen
        self.__spaceCenter = spaceCenter
        self.__ligand = None
        self.__pose_ligand = poseLigand
        self.__typeLS = typeLS
        self.__angleList = angleList


    # -------------------------- Init Process ---------------------------
    def initLocalProcess(self, selectedPop, ligand, alphaGenerations = 1):
        cont = 0
        localT = self.__tempAlpha
        self.__ligand = ligand
        T = self.__temp
        oldSol = copy.deepcopy(selectedPop)
        oldSol.score = self.getScoring(oldSol)
        while( T > self.__tempMin):
            j = 1
            while(j <= self.__innerIteration):
                # Local Search with space reduction by iterations
                if self.__typeLS == 1:
                    alphaRed = T / float(self.__temp)
                    newSol = self.getNeighbor(oldSol, alphaRed)
                # Local Search with space reduction by generations
                elif self.__typeLS == 2:
                    newSol = self.getNeighbor(oldSol, alphaGenerations)
                # Local Search with space reduction by iterations and generations
                elif self.__typeLS == 3:
                    alphaTemp = T / float(self.__temp)
                    alphaRed = alphaGenerations * alphaTemp
                    newSol = self.getNeighbor(oldSol, alphaRed)
                else:
                    newSol = self.getNeighbor(oldSol)
                if(self.acceptance(oldSol, newSol, T) > random.random()):
                    oldSol = copy.deepcopy(newSol)
                j += 1
            T *= self.__tempAlpha
            cont= cont +1
        selectedPop = copy.deepcopy(oldSol)
        return selectedPop
    # --------------------------- Acceptance ----------------------------
    def acceptance(self, oldSol, newSol, temp):
        scoreOld = oldSol.score
        scoreNew = newSol.score
        if(scoreNew <= scoreOld):
            result = 1
        else:
            delthaE = abs(scoreOld - scoreNew)
            k = (1 / (1 + (delthaE / temp)))
            result = math.exp( - delthaE / (temp * k) )
        return result


    # ----------------------- Gets New Neighbor -------------------------
    def getNeighbor(self, oldSol, alpha = 1):
        solution = copy.deepcopy(oldSol)
        solution = self.mutation(solution, alpha)
        solution.score = self.getScoring(solution)
        return solution
    # ----------------------- Search angle list -------------------------

    def searchAngle(self, aList, value, p):
        tempL = list()
	i = 0
	for par in aList:
		if int(par[p]) == value:
			tempL.append(i)
		i = i+1
	return tempL

    #----------------------------- Mutation -----------------------------
    def mutation(self, selectedPop, alpha = 1):
        rotation_bonds = len(selectedPop.rotateAtoms)
        if(rotation_bonds > 0):
            mutPos = random.randint(1, (6 + rotation_bonds))
        else:
             mutPos = random.randint(1, 6)
        if mutPos == 1:
            distright = self.__spaceLen - selectedPop.x
            distleft = -self.__spaceLen + selectedPop.x
            selectedPop.x = random.uniform(distleft*alpha - selectedPop.x, selectedPop.x + distright*alpha)
        elif mutPos == 2:
            distright = self.__spaceLen - selectedPop.y
            distleft = -self.__spaceLen + selectedPop.y
            selectedPop.y = random.uniform(distleft*alpha - selectedPop.y, selectedPop.y + distright*alpha)
        elif mutPos == 3:
            distright = self.__spaceLen - selectedPop.z
            distleft = -self.__spaceLen + selectedPop.z
            selectedPop.z = random.uniform(distleft*alpha - selectedPop.z, selectedPop.z + distright*alpha)
        elif mutPos == 4:
            stheta = round(selectedPop.sph_theta / math.pi * 100.0)
            distright = round((200 - stheta)*alpha)
            distleft = round(stheta*alpha)
            selectedPop.sph_theta = math.pi * random.randint(-distleft + stheta, stheta + distright) / 100.0
        elif mutPos == 5:
            sphi = round(selectedPop.sph_phi / math.pi * 100.0)
            distright = round((100 - sphi)*alpha)
            distleft = round(sphi*alpha)
            selectedPop.sph_phi = math.pi * random.randint(-distleft + sphi, sphi + distright) / 100.0
        elif mutPos == 6:
            theta = round(selectedPop.theta / math.pi * 100.0)
            distright = round((200 - theta)*alpha)
            distleft = round(theta*alpha)
            selectedPop.theta = math.pi * random.randint(-distleft + theta, theta + distright) / 100.0
        elif mutPos > 6:
            angle = round(selectedPop.rotateAtoms[mutPos - 7] / math.pi * 100.0)
            distright = round((200 - angle)*alpha)
            distleft = round(angle*alpha)
            selectedPop.rotateAtoms[mutPos - 7] = math.pi * random.randint(-distleft + angle, angle + distright) / 100.0
        return selectedPop

    '''
    #----------------------- Mutation with space reduction in each iteration-----------------------------
    def mutationReduce(self, selectedPop, alpha):
        rotation_bonds = len(selectedPop.rotateAtoms)
        if(rotation_bonds > 0):
            mutPos = random.randint(1, (6 + rotation_bonds))
        else:
             mutPos = random.randint(1, 6)
        if mutPos == 1:
            selectedPop.x = random.uniform(-self.__spaceLen*alpha, self.__spaceLen*alpha)
        elif mutPos == 2:
            selectedPop.y = random.uniform(-self.__spaceLen*alpha, self.__spaceLen*alpha)
        elif mutPos == 3:
            selectedPop.z = random.uniform(-self.__spaceLen*alpha, self.__spaceLen*alpha)
        elif mutPos == 4:
            selectedPop.sph_theta = math.pi * random.randint(0, 200*alpha) / 100.0
        elif mutPos == 5:
            selectedPop.sph_phi = math.pi * random.randint(0, 100*alpha) / 100.0
        elif mutPos == 6:
            selectedPop.theta = math.pi * random.randint(0, 200*alpha) / 100.0
        elif mutPos > 6:
            selectedPop.rotateAtoms[mutPos - 7] = math.pi * random.randint(0, 200*alpha) / 100.0

        return selectedPop
    '''
    # -------------------------- Gets Scoring ---------------------------
    def getScoring(self, solution):
        auxLigand = copy.deepcopy(self.__ligand)
        auxLigand.traslateToPoint([(self.__spaceCenter[0] + solution.x),
                                   (self.__spaceCenter[1] + solution.y),
                                   (self.__spaceCenter[2] + solution.z)])
        auxLigand.rotateByVector(Methods.spherePoint(1, solution.sph_theta, solution.sph_phi), solution.theta)
        for r in range(len(self.__ligand.rotationsPoints)):
            auxLigand.rotateSegment(r, solution.rotateAtoms[r])

        'Score'
        auxPose = rosetta.core.pose.Pose()
        auxPose.assign(self.__pose_ligand)
        for atom in xrange(auxPose.residue(auxPose.total_residue()).natoms()):
            atomV = rosetta.numeric.xyzVector_Real()
            atomV.x = round(auxLigand.x[atom],3)
            atomV.y = round(auxLigand.y[atom],3)
            atomV.z = round(auxLigand.z[atom],3)
            auxPose.residue(auxPose.total_residue()).set_xyz(auxLigand.atom[atom], atomV)
        result = self.__scoreFx.generateScoringByPose(auxPose, auxLigand, False)

        return result
