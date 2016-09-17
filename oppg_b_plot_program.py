#Program som plotter resultatene skrevet til u-fil

from scitools.std import *
import matplotlib.pyplot as plt

infile = open("../u_file.txt", "r")

u = []
for line in infile:
    words = line.split()
    u.append(float(words[0]))


infile.close()

N = len(u)-2
h = 1.0/(N+1)
x = zeros(N+2)

for i in range(N+2):
    x[i] = h*i



def u_exact(x):
    """
    Den eksakte losningen av problemet
    """

    return 1 - (1-exp(-10))*x-exp(-10*x)

axis_font={"size":"20"}
plt.plot(x, u)
plt.plot(x, u_exact(x))
plt.legend(("numeric", "analytic"))
plt.xlabel("x", **axis_font)
plt.ylabel("Value of u", **axis_font)
plt.title("Numeric vs Analytic Solution,  GeneralSolver, n=1000", **axis_font)
plt.show()

