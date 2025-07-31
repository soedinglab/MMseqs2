import numpy as np
import matplotlib.pyplot as plt

#I_min = 512 + 256
I_min = 512
I_max = 1024
#I_max = 1024 + 256
L = 16
K_min = I_max
K_max = 1 << 20

overshoots = []
I_list = []

for K in range(K_min, K_max + 1, L):
    min_overshoot = 10 * I_max # large number
    best_I = -1

    for I in range(I_min, I_max + 1, L):
        ceil_K = ((K + I - 1) // I) * I
        overshoot = ceil_K - K

        if overshoot <= min_overshoot:
            min_overshoot = overshoot
            best_I = I

    overshoots.append(min_overshoot)
    I_list.append(best_I)

    if K % (K_max // 8) == 0:
        print("K: %d" % K);

overshoots = np.array(overshoots)
I_list = np.array(I_list)
print("Average overshoot: %f" % np.mean(overshoots))
print("Max overshoot: %d" % np.max(overshoots))
print("Average I: %f" % np.mean(I_list))

plt.subplot(121)
plt.hist(overshoots, bins = "auto")
plt.title("Overshoots")

plt.subplot(122)
plt.hist(I_list, bins = "auto")
plt.title("I")

plt.show()
