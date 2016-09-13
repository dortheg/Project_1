#Program som plotter 

from scitools.std import *

infile = open("u_file.txt", "r")

u = []
for line in infile:
    words = line.split()
    u.append(float(words[0]))


infile.close()

N = len(u)
h = 1.0/(N+1)
x = zeros(N)

for i in range(N):
    x[i] = h*i



def u_exact(x):
    """
    Den eksakte losningen av problemet
    """

    return 1 - (1-exp(-10))*x-exp(-10*x)

plot(x,u, legend="numeric")
hold("on")
plot(x, u_exact(x), legend="analytic")
hold("off")
