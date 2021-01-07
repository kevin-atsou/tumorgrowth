from fipy import numerix

# Axi-symmetric form of sigma
def AxiSymmetricSigmaF(r, Rs):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz**2))/(numerix.sqrt(dz**2)*numerix.sqrt(2*numerix.pi)))
    return 10*numerix.exp(-((r)**2)/(2*Rs))/(numerix.sqrt(2*numerix.pi*Rs))

def SigmaF1D(x, Rs):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz**2))/(numerix.sqrt(dz**2)*numerix.sqrt(2*numerix.pi)))
    # return 0.00001*numerix.exp(-((x)**2)/(2*Rs))/(numerix.sqrt(2*numerix.pi*Rs))
    return 10.*numerix.exp(-((x)**2)/(2*Rs))/(numerix.sqrt(2*numerix.pi*Rs))

def SigmaF2D(amp, x,y, Rs):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz**2))/(numerix.sqrt(dz**2)*numerix.sqrt(2*numerix.pi)))
    return amp*numerix.exp(-(x**2 + y**2)/(2*Rs))/(numerix.sqrt(2*numerix.pi*Rs))

def TwoSigmaF2D(x,y, Rs):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz**2))/(numerix.sqrt(dz**2)*numerix.sqrt(2*numerix.pi)))
    return ((200e-3+0.)*numerix.exp(-((x+0.6)**2 + (y-0.3)**2)/(2*Rs))
            +(300e-3+0.)*numerix.exp(-((x-0.6)**2 + (y-0.2)**2)/(2*Rs)))/(numerix.sqrt(2*numerix.pi*Rs))