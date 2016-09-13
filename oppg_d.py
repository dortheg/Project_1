from scitools.std import *
from os import system

N = [10, 100, 1000, 10000]

def u_exact(x):
    """
    Den eksakte losningen av problemet
    """

    return 1 - (1-exp(-10))*x-exp(-10*x)

eps_max = []

for n in N: 
    systemStr = "./Project_1_c %d" % (n)
    system(systemStr)


    x = zeros(n)
    h = 1.0/(n+1)

    for i in range(n):
        x[i] = h*i
    
    u_exact_array = u_exact(x)

    infile = open("u_file.txt", "r")

    u = []
    for line in infile:
        words = line.split()
        u.append(float(words[0]))
    infile.close()
    u = array(u)

    eps = log (abs((u[1:n-1] - u_exact_array[1:n-1])/(u_exact_array[1:n-1])))
    eps_max.append(max(eps))

print eps_max


    
