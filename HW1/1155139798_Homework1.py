#For HW1 - program to integrate using trapezoid and monte Carlo methods

import random, math
import numpy as np
import matplotlib.pyplot as plt


#Analytical result of the integration
trueValue = 0.4


def funcForIntegration(x):
    return x**4


def monteCarloInt(func,M,start,stop):
    '''

        Arguments:

        func -- The function corresponding to the integral being calculated

        M -- The number of sample points

        start -- Lower bound for integration

        stop -- Upper bound for integration
    '''
    s = 0
    ds2 = 0

    if start > stop:
        start, stop = stop, start

    for i in range(0,M):

        x_rand = random.uniform(start,stop)
        s += func(x_rand)
        ds2 += func(x_rand)**2


    s *= (stop - start)/M
    ds2 /= M

    ds2 = math.sqrt(math.fabs((ds2 - s**2))/M)
    return s,ds2


def trapezoidMethod(func,M,start,stop):
    '''

        Arguments:

        func -- The function corresponding to the integral being calculated

        M -- The number of slices

        start -- Lower bound for integration

        stop -- Upper bound for integration
    '''
    if start > stop:
        start, stop = stop, start


    h = (stop - start)/M
    s = 0
    ds2 = 0
    for i in range(0,M):

        s += func(start+h*i) + func(start + h*(i+1))


    s *= h/2

    ds2 = math.fabs(trueValue - s)/trueValue


    return s, ds2


#Number of slices/sample points for calculation
M = 10**6
answerMC = monteCarloInt(funcForIntegration,M,-1,1)
answerTM = trapezoidMethod(funcForIntegration,M,-1,1)


#Create a list of powers of 10
M_values = [10**x for x in range(1,7)]
#List comprehension to give the variance at various M for both methods
varianceMC = [monteCarloInt(funcForIntegration,x,-1,1)[1] for x in M_values]
varianceTM = [trapezoidMethod(funcForIntegration,x,-1,1)[1] for x in M_values]






# In[]:

#Using matplotlib to plot graphs of the data
plt.figure(figsize=(8,8))

plt.subplot(211)
plt.title("Trapezoid Method")
plt.ylabel("Variation")
plt.xlabel("Number of slices")
plt.plot(M_values,varianceTM)
plt.yscale("log")
plt.xscale("log")
plt.grid(True)

plt.subplot(212)
plt.title("Monte Carlo Method")

plt.ylabel("Variation")
plt.xlabel("Number of random data points")
plt.plot(M_values,varianceMC)
plt.yscale("log")
plt.xscale("log")
plt.grid(True)

plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.50,
                    wspace=0.5)


plt.show()






# In[]:


print("The squared deviations for Trapezoidal rule for various M are given by (M,variance):")
print([(x,y) for x,y in zip(M_values,varianceTM)])
print("\n")

print("The squared deviations for Trapezoidal rule for various M are given by:")
print([(x,y) for x,y in zip(M_values,varianceMC)])
print("\n")

print("The variance multiplied by M^2 give a constant value (except when the variance is very low due to floating point errors) as shown below:")
print(np.array(varianceTM)*np.power(np.array(M_values),2))
print("\n")

print("The variance multiplied by M^2 give a roughly constant value as shown below:")
print(np.array(varianceMC)*np.sqrt(np.array(M_values)))
print("\n")



# In[]:

print("The approximation by trapezoidal rule for {0} slices is {1:}".format(M,answerTM[0]))
print("The approximation by trapezoidal rule for {0} slices is {1:}".format(M,answerMC[0]))
