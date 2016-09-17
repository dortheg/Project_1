from scitools.std import *
import matplotlib.pyplot as plt
from os import system

N = [10, 100, 1000, 10000, 100000, 10**6, 10**7]

def u_exact(x):
    """
    Den eksakte losningen av problemet
    """

    return 1 - (1-exp(-10))*x-exp(-10*x)

eps_max = []

for n in N: 
    systemStr = "../build-Project_1-Desktop_Qt_5_7_0_GCC_64bit-Release/Project_1 %d" % (n)
    system(systemStr)


    x = zeros(n+2)
    h = 1.0/(n+1)

    for i in range(n+2):
        x[i] = h*i
    
    u_exact_array = u_exact(x)

    infile = open("../u_file.txt", "r")

    u = []
    for line in infile:
        words = line.split()
        u.append(float(words[0]))
    infile.close()
    u = array(u)
    
    #kan bare sammenlikne med verdiene u[1:n], siden endepunktene er 0
    eps =  log10 (abs((u[1:n] - u_exact_array[1:n])/(u_exact_array[1:n])))
    eps_max.append(max(eps))

outfile = open("../eps_file.txt", "w")
for epsilon in eps_max:
    outfile.write("%.5f \n" % epsilon)
outfile.close()

axis_font={"size":"20"}
plt.plot(log10(N), eps_max)
plt.xlabel("Number of grid points n, log_10 scale", **axis_font)
plt.ylabel("Relative error, log_10 scale", **axis_font)
plt.title("Relative error, Numeric vs Analytic Solution", **axis_font)
plt.show()
    
