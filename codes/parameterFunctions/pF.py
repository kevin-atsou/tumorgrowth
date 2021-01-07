from fipy import numerix

# Axi-symmetric form of sigma
def AxiSymmetricGaussianConversionRate(r, mean, Rp):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz**2))/(numerix.sqrt(dz**2)*numerix.sqrt(2*numerix.pi)))
    return numerix.exp(-(r-mean)**2/(2*Rp))/(numerix.sqrt(2*numerix.pi*Rp))

def GaussianConversionRate(r, mean, Rp):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz**2))/(numerix.sqrt(dz**2)*numerix.sqrt(2*numerix.pi)))
    return 0.3*numerix.exp(-(r-mean)**2/(2*Rp))/(numerix.sqrt(2*numerix.pi*Rp))

def GaussianConversionRate2D(amp, x,y, mean, Rp):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz**2))/(numerix.sqrt(dz**2)*numerix.sqrt(2*numerix.pi)))
    return amp*numerix.exp(-((x-mean)**2 + (y-mean)**2)/(2*Rp))/(numerix.sqrt(2*numerix.pi*Rp))