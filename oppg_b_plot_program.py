#Program som plotter 

from scitools.std import *
import matplotlib.pyplot as plt

infile = open("u_file.txt", "r")

u = []
for line in infile:
    words = line.split()
    u.append(float(words[0]))


infile.close()

N = len(u)-2
h = 1.0/(N+1)
x = zeros(N+2)

#Denne x-en er en annen i programmet, siden denne skal sendes til u_exact
for i in range(N+2):
    x[i] = h*i



def u_exact(x):
    """
    Den eksakte losningen av problemet
    """

    return 1 - (1-exp(-10))*x-exp(-10*x)

"""
plot(x,u, legend="numeric")
hold("on")
plot(x, u_exact(x), legend="analytic")
hold("off")
"""

axis_font={"size":"20"}
plt.plot(x, u)
plt.plot(x, u_exact(x))
plt.legend(("numeric", "analytic"))
plt.xlabel("x", **axis_font)
plt.ylabel("Value of u", **axis_font)
plt.title("Numeric vs Analytic Solution,  GeneralSolver, n=10", **axis_font)
plt.show()

