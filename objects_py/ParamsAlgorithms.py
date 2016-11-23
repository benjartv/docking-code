

# ------------------------ Params Algorithms ------------------------
class ParamsAlgorithms:

    def __init__(self):
        self.searchSpaceSize = 0
        self.searchCenterPoint = []
        self.newCenterPoint = []
        self.moleculeName = ""
        self.testPdb = ""
        self.isLocalSearch = False
        self.rotatesAtoms = []
        self.populationSize = 0
        self.generationNumber = 0
        self.populationCast = []
        self.algorithmType = 0
        self.folderName = ""
        self.timeStop = 0
        self.loopReset = 0
        self.mutProbability = 0
        self.typeLS = 0
        self.typeCO = 0
        self.distanceCriteria = 0
        self.angleList = []

        'Local Search'
        self.initTemperature = 0
        self.minTemperature = 0
        self.alphaTemperature = 0
        self.innerGenerations = 0
        self.angleList = []
