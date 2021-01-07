from fipy import numerix

# Axi-symmetric form of sigma
def gaussianDivisionRate(z,mean,sigma):
    #Test : plt.plot(r,10*numerix.exp(-r**2/(2*dz))/numerix.sqrt(2*numerix.pi*2*dz))
    return numerix.exp(-(z-mean)**2/(2*sigma**2))/(sigma*numerix.sqrt(2*numerix.pi))