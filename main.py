import numpy as np
import matplotlib.pyplot as plt
from zkrp import frechet_mean
from zkrp import hausdorff_distance
from zkrp import hausdorff_distance2
from zkrp import frechet_variance
from zkrp import frechet_covariance
from zkrp import frechet_covariance2
from zkrp import dist
from zkrp import show

np.random.seed(0)
n = 100

central = np.append(np.random.normal(0, 1, round(n / 4)), np.random.normal(4, 1, round(n * 3 / 4)))
# central = np.random.normal(0, 1, n)

mydata = np.zeros((n, 2))
for i in range(n):
    interval_len = abs(np.random.normal(0, 1, 1))
    mydata[i, 0] = central[i] - interval_len / 2
    mydata[i, 1] = central[i] + interval_len / 2

# show(mydata, './fig/interval data bimodal.jpg')

result = frechet_mean(mydata, method='hausdorff')
result2 = frechet_mean(mydata, method='middle')
# result_list = []
# length = [0.2, 0.5, 1, 1.5, 2]
# color = ['red', 'green', 'blue', 'orange', 'black']
# for i in range(len(length)):
#     result_list.append(frechet_mean(mydata, method='hausdorff', inter_len=length[i]))


# show(mydata, './fig/interval data with mean.jpg', result)
# show(mydata, './fig/interval data bimodal with mean3.jpg', result)
# for i in range(len(length)):
#     dsum = result_list[i]['dsum']
#     x = range(len(dsum))
#     plt.plot(x, dsum, alpha=1, color=color[i], label='interval length: %.1f' % length[i])

dsum = result['dsum']
x = range(len(dsum))
dsum2 = result2['dsum']

plt.plot(x, dsum, alpha=1, color='red', label='hausdorff')
plt.plot(x, dsum2, alpha=1, color='blue', label='middle')
# plt.title('plot of normal distribution middle')  # 折线图标题
plt.title('plot of bimodal distribution with diff method')  # 折线图标题

plt.xlabel('X')  # x轴标题
plt.ylabel('Sum d')  # y轴标题
plt.grid(ls='-.')  # 绘制背景线
plt.legend(loc='best')
plt.tight_layout()
plt.show()

# plt.savefig('./fig/normal distribution.jpg')
plt.savefig('./fig/bimodal distribution with diff method.jpg')

plt.close()

# =====================================
np.random.seed(0)
n = 100

centralx = np.linspace(1, 5, n)
centraly = np.linspace(5, 1, n)
centralz = np.random.normal(3, 1, n)

mydatax = np.zeros((n, 2))
interval_lenx = abs(np.random.normal(0, 1, n))
mydatax[:, 0] = centralx - interval_lenx / 2
mydatax[:, 1] = centralx + interval_lenx / 2

mydatay = np.zeros((n, 2))
interval_leny = abs(np.random.normal(0, 1, n))
mydatay[:, 0] = centraly - interval_lenx / 2
mydatay[:, 1] = centraly + interval_lenx / 2

for i in range(n):
    print(hausdorff_distance(mydatax[i,], mydatay[i,]) == hausdorff_distance2(mydatax[i,], mydatay[i,]))

np.random.seed(0)
n = 100

centralx = np.linspace(1,5,n)
centraly = np.linspace(5,1,n)
centralz = np.random.normal(3, 1, n)

mydatax = np.zeros((n, 2))
interval_lenx = abs(np.random.normal(0, 1, n))
mydatax[:, 0] = centralx - interval_lenx / 2
mydatax[:, 1] = centralx + interval_lenx / 2

data = mydatax
i = 491
lower = min(data[:, 0])
upper = max(data[:, 1])
inter_len = np.mean(data[:, 1] - data[:, 0])
d_l = []
n = np.shape(data)[0]


d = 0
interval = [lower + (upper - lower) * i / 1000 - inter_len / 2,
            lower + (upper - lower) * i / 1000 + inter_len / 2]
for j in range(n):
    d_l.append(hausdorff_distance(interval, data[j, :]) * hausdorff_distance(interval, data[j, :]))
    d += hausdorff_distance(interval, data[j, :]) * hausdorff_distance(interval, data[j, :])

dsum.append(d)
min_d = min(dsum)
min_index = dsum.index(min_d)
min_interval = [lower + (upper - lower) * min_index / 1000 - 0.5,
                lower + (upper - lower) * min_index / 1000 + 0.5]