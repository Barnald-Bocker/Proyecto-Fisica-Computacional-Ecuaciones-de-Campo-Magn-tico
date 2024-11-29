#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

file = pd.read_csv('capacitor.csv',delimiter=';', header=None)

data = file.to_numpy(dtype=float)
plt.imshow(data)
plt.gray()
plt.show()
