import numpy as np
import matplotlib.pyplot as plt

def find_up_down(x, y, p):
        up = np.min(x[x>p])
        upy = np.min(y[x>p])
        down = np.max(x[x<p])
        downy = np.max(y[x<p])
        px = p
        m = (upy - downy) / (up - down)
        py = (m * (px - down)) + downy
        return up, upy, down, downy, px, py


x = np.linspace(0, 40, 40)


x1 =  np.linspace(0, 10, 10)
x2 =  np.linspace(11, 20, 10)
x3 =  np.linspace(21, 30, 10)
x4 =  np.linspace(31, 40, 10)

y1 = x1
y2 = (4/3)* x2 + 3
y3 = np.pi * 0.5 * x3 + 0.5
y4 = x4 + 6

y = np.hstack([y1, y2, y3, y4])


print(len(x))
print(len(y))
plt.subplot()
plt.scatter(x, y)
plt.plot(x,y)
p = 13.1
up, upy, down, downy, px, py = find_up_down(x, y,p)
plt.scatter([up, down], [upy, downy])
plt.scatter(px, py)


plt.grid(color='white', ls='solid', linewidth = 0.3)

plt.show()



