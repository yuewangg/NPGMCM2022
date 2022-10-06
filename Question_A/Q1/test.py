import numpy as np
import os

pwd = os.getcwd()
print(pwd)
loadData = np.load('../data_q1.npy')

print("----type----")
print(type(loadData))
print("----shape----")
print(loadData.shape)