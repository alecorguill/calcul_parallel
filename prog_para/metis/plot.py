import matplotlib.pyplot as plt

nb_procs, times = [], []
with open("time.plt") as fp:
    for line in fp:
        tab = line.split(" ")
        nb_procs.append(int(tab[0]))
        times.append(float(tab[1]))

times = [times[0]/val for val in times]
plt.plot(nb_procs, times)
plt.show()
