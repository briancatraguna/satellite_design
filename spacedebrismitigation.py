import matplotlib.pyplot as plt

Vinit = 7672.594
Vfin = 7769.480

def calc_height(v):
    height = (3.986*10**14)/(v**2)
    return height

def convertmtokm(m):
    return m/1000

def convertstodays(s):
    return (((s/60)/60)/24)

v = [Vinit]
t=[0]
i=0
acc=float(1.259*10**(-4))
while v[i]<=Vfin:
    t.append(t[i]+1)
    v.append(v[i]+acc)
    i+=1

height = []
for i,vel in enumerate(v):
    height.append(calc_height(vel))
    height[i] = convertmtokm(height[i])-6371
    t[i] = convertstodays(t[i])

plt.title("Altitude over Time")
plt.plot(t,height)
plt.xlabel("Time (days)")
plt.ylabel("Height (km)")
plt.show()

print("Time to re-enter:",t[-1])