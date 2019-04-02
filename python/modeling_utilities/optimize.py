# minimize a function computed by another code

import subprocess
import numpy as np
from scipy.optimize import minimize

def func(x):
    s = subprocess.check_output(["python", "/Users/cerusu/GITHUB/zMstarPDF/python/modeling_utilities/func.py", str(x[0]), str(x[1])])
    print str(x[0]), str(x[1])
    return float(s)

x0 = np.array([0.1,0.1])
res = minimize(func, x0, method='nelder-mead', options={'xtol': 1e-8, 'disp': True})
print(res.x)
