from fipy import numerix
import scipy.sparse as sp;
import  numpy as np

def gaussian(x, mu, sig):
    return numerix.exp(-numerix.power(x - mu, 2.) / (2 * numerix.power(sig, 2.)))

def loopGaussianKernel(z, dz):
    n = z.shape[0]
    m = z.shape[0]
    kern = sp.lil_matrix((n,m), dtype=np.float64)
    # kern = numerix.zeros((n,m))
    # eps = 1e-10
    eps = dz**2
    digits = len(str(n-1))
    delete = "\b" * (digits + 1)
    for i in range(n):
        #if i>0:
            # for j in range(m):
               # if j>0:
               #      kern[i, j] = numerix.exp((-(2*z[i] - z[j])**2)/(2*eps))
        kern[i] = numerix.exp((-(2*z[i] - z)**2)/(2*eps))
                    # kern = kern + sp.coo_matrix((np.array([numerix.exp((-(2*z[i] - z[j])**2)/(2*eps))]), (np.array([i]), np.array([j]))), shape=(n, m))
                    # kern.tocsc()
                    #/numerix.sqrt(2*numerix.pi*eps)
            # print "{0}{1:{2}}".format(delete, i, digits),

        # level = int((np.float(i)/np.float(n))*100)
        # print "{0}%\r".format((np.float(i)/np.float(n)*100))
        # print "#"*(level+1),
    return kern/numerix.sqrt(2*numerix.pi*eps)


def gaussianKernel(z, dz):
    Z = np.ones((z.shape[0], 1)) * np.array([z])
    eps = dz**2
    kern = np.exp(-(2 * Z.T - Z) * (2 * Z.T - Z) / (2 * eps))
    return kern/numerix.sqrt(2*numerix.pi*eps)


def UniformKernel(z):
    n = z.shape[0]
    m = z.shape[0]
    kern = sp.lil_matrix((n,m), dtype=np.float64)
    temp = numerix.zeros(n)
    digits = len(str(n-1))
    delete = "\b" * (digits + 1)
    for i in range(n):
        indices = np.where(z >= z[i])
        temp[indices] = 1./z[z >= z[i]]
        kern[i] = temp
        temp = numerix.zeros(n)
        level = int((np.float(i)/np.float(n))*100)
        # print "{0}%\r".format((np.float(i)/np.float(n)*100))
    return kern
