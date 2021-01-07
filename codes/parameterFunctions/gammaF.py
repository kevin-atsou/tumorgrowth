from fipy import numerix

# Axi-symmetric form of sigma
def AxiSymmetricGaussianGammaF(r, Rs):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz**2))/(numerix.sqrt(dz**2)*numerix.sqrt(2*numerix.pi)))
    return 0.03*numerix.exp(-(r**2)/(2*Rs))/(numerix.sqrt(2*numerix.pi*Rs))

def gaussianGammaF(x,Rs):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz**2))/(numerix.sqrt(dz**2)*numerix.sqrt(2*numerix.pi)))
    return 0.3*numerix.exp(-(x**2)/(2*Rs))/(numerix.sqrt(2*numerix.pi*Rs))

def gaussianGammaF2D(amp,x,y,Rs):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz**2))/(numerix.sqrt(dz**2)*numerix.sqrt(2*numerix.pi)))
    return amp*numerix.exp(-(x**2 + y**2)/(2*Rs))/(numerix.sqrt(2*numerix.pi*Rs))