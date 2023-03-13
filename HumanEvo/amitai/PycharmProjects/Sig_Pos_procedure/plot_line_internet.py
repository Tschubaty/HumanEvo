import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy import stats

interv_names = []

X_interv = []
interv_num_samples = []
INTERVAL = 1 #
y = np.array([0.92271, 1.,      1.,      1.,      1. ,     1. ,     0.83354, 1. ,     1.
, 1.,      1. ,     0.85678 ,1.  ,    1. ,     1.    ,  1.    ,  1.  ,    0.98873,
 0.93889, 0.78555 ,0.64624, 0.84454, 0.9991 , 1.   ,   1.   ,   1.  ,    0.90566,
 0.99679, 1. ,     1.   ,   1.   ,   0.66173, 1. ,     1.   ,   0.85443, 0.88852,
 1. ,     1.   ,   1.   ,   1.    ,  0.95761 ,1.  ,    0.89348 ,1.,      0.81574,
 0.54157 ,0.67293 ,0.78095, 1,      1.   ,   0.9961 , 1.  ,    0.81742 ,1.,
 0.75575, 1.,      1.   ,   0.63144, 0.94447 ,0.72199, 0.89982 ])

X = np.array(  [1, 1, 1.2, 1.6, 1.8, 2.2, 2.5, 2.8, 3, 3.1, 3.3, 3.4, 3.5, 3.8, 3.8, 4, 4.2, 4.3, 4.4, 4.5, 4.5, 4.6, 4.7, 4.7,
     4.9, 5.1, 5.1, 5.1, 5.2, 5.4, 5.4, 5.8, 6.1, 6.1, 6.1, 7, 7, 7, 7.3, 7.3, 7.5, 7.6, 7.7, 7.8, 7.9, 7.9, 8, 8, 8, 8,
     8.1, 8.1, 8.1, 8.4, 8.5, 8.9, 9, 9.2, 10, 10.8, 12])



def create_intervals(times):
    cur_interv = 1
    num_samples_cur_interval = 0
    for i, time in enumerate(times):
        if((time< cur_interv + INTERVAL) and (time>= cur_interv)):

            num_samples_cur_interval +=1
        else:
            if(num_samples_cur_interval != 0):
                interv_names.append(str(cur_interv) + '-' + str(cur_interv + INTERVAL) + 'K_'+str(num_samples_cur_interval) + '_samples')
                X_interv.append(np.mean(times[i - num_samples_cur_interval: i]))
                interv_num_samples.append(num_samples_cur_interval)
                cur_interv = int(np.min(times[i:]))
                num_samples_cur_interval = 1

        if (i == len(times) - 1):
            interv_names.append(str(cur_interv) + '-' + str(cur_interv + INTERVAL) + 'K_' + str(num_samples_cur_interval) + '_samples')
            X_interv.append(np.mean(time))
            interv_num_samples.append(1)


create_intervals(X)
y_interv= []
print(interv_num_samples)
print(X_interv)

cur_loc = 0
for elem in interv_num_samples:
    y_interv.append(np.mean(y[cur_loc:cur_loc + elem]))
    cur_loc += elem
print(y_interv)

y_interv = np.array(y_interv)
X_interv = np.array(X_interv)




def draw_line(X,y):
    res_arr  = stats.linregress(X, y)


    line_y= res_arr[1] + res_arr[0]*X
    # plot the results
    print(res_arr)
    plt.figure(figsize=(4, 3))
    ax = plt.axes()
    ax.scatter(X, y)
    ax.plot(X, line_y, color='red')

    ax.set_xlabel('Methylation Level')
    ax.set_ylabel('y')

    ax.axis('tight')


    plt.show()



print(X)
print(y)
print(X_interv)
print(y_interv)
draw_line(X_interv,y_interv)