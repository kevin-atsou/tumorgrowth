from fipy import numerix


def GaussianKernel(strength, x, Ra):
    return strength*numerix.exp(-(x**2)/(2*Ra))/(numerix.sqrt(2*numerix.pi*Ra))


def GaussianKernel2D(strength, x, y, Ra):
    return strength*numerix.exp(-((x**2)+(y**2))/(2*Ra))/(numerix.sqrt(2*numerix.pi*Ra))