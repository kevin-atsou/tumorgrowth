class TestFeature:
    def __init__(self, a, V, delta, pf, S, chi, D, sF, gf, k, tested_param, init_val, alpha, beta,):
        self.a = a
        self.V = V
        self.delta = delta
        self.pf = pf
        self.S = S
        self.chi = chi
        self.D = D
        self.sF = sF
        self.gf = gf
        self.k = k
        self.tested_param = tested_param
        self.init_val = init_val
        self.alpha = alpha
        self.beta = beta

    def setA(self, value):
        self.a = value

    def getA(self):
        return self.a

    def setV(self, value):
        self.V = value

    def getV(self):
        return self.V

    def setDelta(self, value):
        self.delta = value

    def getDelta(self):
        return self.delta

    def setPf(self, value):
        self.pf = value

    def getPf(self):
        return self.pf

    def setS(self, value):
        self.S = value

    def getS(self):
        return self.S

    def setChi(self, value):
        self.Chi = value

    def getChi(self):
        return self.Chi

    def setD(self, value):
        self.D = value

    def getD(self):
        return self.D

    def setSf(self, value):
        self.Sf = value

    def getSf(self):
        return self.Sf

    def setGf(self, value):
        self.gf = value

    def getGf(self):
        return self.gf

    def setK(self, value):
        self.k = value

    def getK(self):
        return self.k

    def setInit_Value(self, value):
        self.init_val = value

    def getInit_Value(self):
        self.init_val

    def setAlpha(self, value):
        self.alpha = value

    def getAlpha(self):
        return self.alpha

    def setBeta(self, value):
        self.beta = value

    def getBeta(self):
        return self.beta

    def setTestedParam(self, value):
        self.tested_param = value

    def getTestedParam(self):
        return self.tested_param

    def setAttribute(self, i, value):
        if i == 0:
            self.a = value
            self.tested_param = 'a'
        if i == 1:
            self.V = value
            self.tested_param = 'V'
        if i == 2:
            self.delta = value
            self.tested_param = 'delta'
        if i == 3:
            self.pf = value
            self.tested_param = 'pf'
        if i == 4:
            self.S = value
            self.tested_param = 'S'
        if i == 5:
            self.chi = value
            self.tested_param = 'chi'
        if i == 6:
            self.D = value
            self.tested_param = 'D'
        if i == 7:
            self.sF = value
            self.tested_param = 'sF'
        if i == 8:
            self.gf = value
            self.tested_param = 'gf'
        if i == 9:
            self.k = value
            self.tested_param = 'd'
        if i == 10:
            self.init_val = value
            self.tested_param = 'init'
        if i == 11:
            self.alpha = value
            self.tested_param = 'alpha'
        if i == 12:
            self.beta = value
            self.tested_param = 'beta'

    def getAttributeByIndex(self, i):
        if i == 0:
            return self.a
        if i == 1:
            return self.V
        if i == 2:
            return self.delta
        if i == 3:
            return self.pf
        if i == 4:
            return self.S
        if i == 5:
            return self.chi
        if i == 6:
            return self.D
        if i == 7:
            return self.sF
        if i == 8:
            return self.gf
        if i == 9:
            return self.k
        if i == 10:
            return self.init_val
        if i == 11:
            return self.alpha
        if i == 12:
            return self.beta

    def getAttributeByName(self, att):
        if att == 'a':
            return self.a
        if att == 'V':
            return self.V
        if att == 'delta':
            return self.delta
        if att == 'pf':
            return self.pf
        if att == 'S':
            return self.S
        if att == 'chi':
            return self.chi
        if att == 'D':
            return self.D
        if att == 'sF':
            return self.sF
        if att == 'gf':
            return self.gf
        if att == 'd':
            return self.k
        if att == 'init':
            return self.init_val
        if att == 'alpha':
            return self.alpha
        if att == 'beta':
            return self.beta

    def setAttributeByName(self, att, value):
        if att == 'a':
            self.a = value
        if att == 'V':
            self.V = value
        if att == 'delta':
            self.delta = value
        if att == 'pf':
            self.pf = value
        if att == 'S':
            self.S = value
        if att == 'chi':
            self.chi = value
        if att == 'D':
            self.D = value
        if att == 'sF':
            self.sF = value
        if att == 'gf':
            self.gf = value
        if att == 'd':
            self.k = value
        if att == 'init':
            self.init_val = value
        if att == 'alpha':
            self.alpha = value
        if att == 'beta':
            self.beta = value


    def getAttributeNameByIndex(self, i):
        if i == 0:
            return 'a'
        if i == 1:
            return 'V'
        if i == 2:
            return 'delta'
        if i == 3:
            return 'pf'
        if i == 4:
            return 'S'
        if i == 5:
            return 'chi'
        if i == 6:
            return 'D'
        if i == 7:
            return 'sF'
        if i == 8:
            return 'gf'
        if i == 9:
            return 'd'
        if i == 10:
            return 'init'
        if i == 11:
            return 'alpha'
        if i == 12:
            return 'beta'

    def getIndexByAttributeName(self, att):

        if att == 'a':
            return 0
        if att == 'V':
            return 1
        if att == 'delta':
            return 2
        if att == 'pf':
            return 3
        if att == 'S':
            return 4
        if att == 'chi':
            return 5
        if att == 'D':
            return 6
        if att == 'sF':
            return 7
        if att == 'gf':
            return 8
        if att == 'd':
            return 9
        if att == 'init':
            return 10
        if att == 'alpha':
            return 11
        if att == 'beta':
            return 12