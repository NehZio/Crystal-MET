import numpy as np
import scipy.optimize as optimize

def func(param,charge,const):
    result = 0
    for i in const:
        result += i[0]*i[1]

    for i in range(len(param)):
        result += param[i]*coef[i]

    return result**2

f = open('toto','r')

charge = []
const = []
coef = []

line = f.readline()

while line:
    ini = 0
    line = float(line.split()[4])
    if line == -2.0 or line == 2.0 or line == 4.0:
        for i in const:
            if i[0] == line:
                i[1] += 1
                ini = 1
        if ini == 0:
                const.append([line,1])
        line = f.readline()
        continue

    for i in range(len(charge)):
        if charge[i] == line:
            coef[i] += 1
            ini = 1
    if ini == 0:
        charge.append(line)
        coef.append(1)

    line = f.readline()
f.close()
results = optimize.minimize(func,charge,args=(coef,const))
#if results.success:
fitted_param = results.x
for i in range(len(charge)):
    print('% 7.5f   =>   % 7.5f  \n'%(charge[i],fitted_param[i]))
#else:
#    print(charge)
#    raise ValueError(results.message)

print('%  7.5f     =>     % 7.5f  \n'%(np.sqrt(func(charge,coef,const)),np.sqrt(func(fitted_param,coef,const))))


