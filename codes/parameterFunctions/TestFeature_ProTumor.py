class TestFeature:
    def __init__(self, a, V, delta, pf, S, chi, D, sF, gf, k, b1, b2, eta, sF2, tau, gammar, kr, theta1, theta2, pfr, kc
                 ,mc, tested_param, init_val, alpha, beta, Sr):
        # Parameters for the effector Cells
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

        # Parameters for the Pro-tumor immune cells
        self.b1 = b1
        self.b2 = b2
        self.eta = eta
        self.sF2 = sF2
        self.tau = tau
        self.gammar = gammar
        self.kr = kr
        self.theta1 = theta1
        self.theta2 = theta2
        self.pfr = pfr
        self.kc = kc
        self.mc = mc

        self.tested_param = tested_param
        self.init_val = init_val
        self.alpha = alpha
        self.beta = beta
        self.Sr = Sr

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

    def setB1(self,value):
        self.b1 = value

    def getB1(self):
        return self.b1

    def setB1(self,value):
        self.b2 = value

    def getB1(self):
        return self.b2

    def setEta(self, value):
        self.eta = value

    def getEta(self):
        return self.eta

    def setSF2(self, value):
        self.sF2 = value

    def getSF2(self):
        return self.sF2

    def setTau(self, value):
        self.tau = value

    def getTau(self):
        return self.tau

    def setGammar(self, value):
        self.gammar = value

    def getGammar(self):
        return self.gammar

    def setKr(self, value):
        self.kr = value

    def getKr(self):
        return self.kr

    def setTheta(self, value):
        self.theta1 = value

    def getTheta(self):
        return self.theta1

    def setTheta(self, value):
        self.theta2 = value

    def getTheta(self):
        return self.theta2

    def setPfr(self, value):
        self.pfr = value

    def getPfr(self):
        return self.pfr

    def setKc(self, value):
        self.kc = value

    def getKc(self):
        return self.kc

    def setMc(self, value):
        self.kc = value

    def getMc(self):
        return self.kc

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

    def setSr(self, value):
        self.S = value

    def getSr(self):
        return self.S

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
            self.b1 = value
            self.tested_param = 'b1'
        if i == 11:
            self.b2 = value
            self.tested_param = 'b2'
        if i == 12:
            self.eta = value
            self.tested_param = 'eta'
        if i == 13:
            self.sF2 = value
            self.tested_param = 'sF2'
        if i == 14:
            self.tau = value
            self.tested_param = 'tau'
        if i == 15:
            self.gammar = value
            self.tested_param = 'gammar'
        if i == 16:
            self.kr = value
            self.tested_param = 'kr'
        if i == 17:
            self.theta1 = value
            self.tested_param = 'theta1'
        if i == 18:
            self.theta2 = value
            self.tested_param = 'theta2'
        if i == 19:
            self.pfr = value
            self.tested_param = 'pfr'
        if i == 20:
            self.kc = value
            self.tested_param = 'kc'
        if i == 21:
            self.mc = value
            self.tested_param = 'mc'
        if i == 22:
            self.init_val = value
            self.tested_param = 'init'
        if i == 23:
            self.alpha = value
            self.tested_param = 'alpha'
        if i == 24:
            self.beta = value
            self.tested_param = 'beta'
        if i == 25:
            self.Sr = value
            self.tested_param = 'Sr'

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
            return self.b1
        if i == 11:
            return self.b2
        if i == 12:
            return self.eta
        if i == 13:
            return self.sF2
        if i == 14:
            return self.tau
        if i == 15:
            return self.gammar
        if i == 16:
            return self.kr
        if i == 17:
            return self.theta1
        if i == 18:
            return self.theta2
        if i == 19:
            return self.pfr
        if i == 20:
            return self.kc
        if i == 21:
            return self.mc
        if i == 22:
            return self.init_val
        if i == 23:
            return self.alpha
        if i == 24:
            return self.beta
        if i == 25:
            return self.S

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
        if att == 'b1':
            return self.b1
        if att == 'b2':
            return self.b2
        if att == 'eta':
            return self.eta
        if att == 'sF2':
            return self.sF2
        if att == 'tau':
            return self.tau
        if att == 'gammar':
            return self.gammar
        if att == 'kr':
            return self.kr
        if att == 'theta1':
            return self.theta1
        if att == 'theta2':
            return self.theta2
        if att == 'pfr':
            return self.pfr
        if att == 'kc':
            return self.kc
        if att == 'mc':
            return self.mc
        if att == 'init':
            return self.init_val
        if att == 'alpha':
            return self.alpha
        if att == 'beta':
            return self.beta
        if att == 'Sr':
            return self.Sr

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
        if att == 'b1':
            self.b1 = value
        if att == 'b2':
            self.b1 = value
        if att == 'eta':
            self.eta = value
        if att == 'sF2':
            self.sF2 = value
        if att == 'tau':
            self.tau = value
        if att == 'gammar':
            self.gammar = value
        if att == 'kr':
            self.kr = value
        if att == 'theta1':
            self.theta1 = value
        if att == 'theta2':
            self.theta2 = value
        if att == 'pfr':
            self.pfr = value
        if att == 'kc':
            self.kc = value
        if att == 'mc':
            self.mc = value
        if att == 'init':
            self.init_val = value
        if att == 'alpha':
            self.alpha = value
        if att == 'beta':
            self.beta = value
        if att == 'Sr':
            self.Sr = value

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
            return 'b1'
        if i == 11:
            return 'b2'
        if i == 12:
            return 'eta'
        if i == 13:
            return 'sF2'
        if i == 14:
            return 'tau'
        if i == 15:
            return 'gammar'
        if i == 16:
            return 'kr'
        if i == 17:
            return 'theta1'
        if i == 18:
            return 'theta2'
        if i == 19:
            return 'pfr'
        if i == 20:
            return 'kc'
        if i == 21:
            return 'mc'
        if i == 22:
            return 'init'
        if i == 23:
            return 'alpha'
        if i == 24:
            return 'beta'
        if i == 25:
            return 'Sr'

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
        if att == 'b1':
            return 10
        if att == 'b2':
            return 11
        if att == 'eta':
            return 12
        if att == 'sF2':
            return 13
        if att == 'tau':
            return 14
        if att == 'gammar':
            return 15
        if att == 'kr':
            return 16
        if att == 'theta1':
            return 17
        if att == 'theta2':
            return 18
        if att == 'pfr':
            return 19
        if att == 'kc':
            return 20
        if att == 'mc':
            return 21
        if att == 'init':
            return 22
        if att == 'alpha':
            return 23
        if att == 'beta':
            return 24
        if att == 'Sr':
            return 25

