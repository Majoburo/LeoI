import numpy as np
import toolbox
import configobj
import configobj
from scipy import interpolate
from matplotlib import pyplot as plt
config = configobj.ConfigObj("config.ini")

data,range= toolbox.getdata(config['vw']['template'])
newrange = np.linspace(range[0],range[-1],num=100000)
newdata = interpolate.InterpolatedUnivariateSpline(range,data)
plt.plot(newrange,newdata(newrange))
plt.show()
