from fipy import numerix

# Axi-symmetric form of sigma
def AxiSymetricGaussianImmuneResponse(r, Ra):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz))/numerix.sqrt(2*numerix.pi*2*dz))
    return 22*numerix.exp(-r**2/(2*Ra))/(numerix.sqrt(2*numerix.pi*Ra))

def GaussianImmuneResponse(strength,x, Ra):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz))/numerix.sqrt(2*numerix.pi*2*dz))
    return strength*numerix.exp(-(x**2)/(2*Ra))/(numerix.sqrt(2*numerix.pi*Ra))

def GaussianImmuneResponse2D(strength,x,y, Ra):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz))/numerix.sqrt(2*numerix.pi*2*dz))
    return strength*numerix.exp(-((x**2)+(y**2))/(2*Ra))/(numerix.sqrt(2*numerix.pi*Ra))

def TwoGaussianImmuneResponse2D(strength,x,y, Ra):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz))/numerix.sqrt(2*numerix.pi*2*dz))
    return ((strength)*numerix.exp(-((x+0.6)**2 + (y-0.3)**2)/(2*Ra))
            +(strength)*numerix.exp(-((x-0.6)**2 + (y-0.2)**2)/(2*Ra)))/(numerix.sqrt(2*numerix.pi*Ra))
