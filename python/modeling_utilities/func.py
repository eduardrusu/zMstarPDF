# outputs the result of a function to stdout
# I use it as the first step to learning how to combine C++ and python code
# run as python func.py 0.1 0.1

import sys
import os
from os import system
import numpy as np

x = float(sys.argv[1])
y = float(sys.argv[2])

print x**2 + y**2
