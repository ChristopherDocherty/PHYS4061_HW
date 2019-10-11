import numpy as np

matA = np.array([[5,2,4,6],[5,8,7,2],[6,4,8,5],[5,4,1,3]])
matB = np.array([[55,23,1,45],[5,65,14,26],[25,14,88,2],[25,25,18,99]])



print("This is the inverse of matrix matA: \n")
print(np.linalg.inv(matA))

print("\n This is the inverse of matrix matB: \n")
print(np.linalg.inv(matB))
