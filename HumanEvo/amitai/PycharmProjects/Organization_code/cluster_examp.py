print("hello world")

import numpy
print("maccabi")

import sys
import csv





import pandas as pd

x = pd.read_csv(r'example.csv')
print(x)

sys.path.insert(0, r'Organization_code')
m = __import__('Arrays_and_dictionaries', globals(), locals(), [])
del sys.path[0]


import matplotlib.pyplot as plt

plt.scatter([1,2,3], [1,2,3])
plt.savefig('graph.png')



print(m.NAN)
