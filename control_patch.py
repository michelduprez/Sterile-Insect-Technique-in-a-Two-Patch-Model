from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from gekko import GEKKO
plt.rcParams['font.size'] = 15



"""
Optimal control of the TIS
model from the paper Bossin-Dumont-Strugarek

"""

# Model arameters
b1 = 8.15
b2 = 8.15
mu11 = 0.035
mu12 = 0.035
mu21 = 0.0015
mu22 = 0.0014
gamma = 1.0
mu_s = 0.233
omega = 1.0


# figure
figure = 4
if figure == 4:
    # migration parameters
    alpha12 = 0.0
    alpha21 = 0.0
if figure == 5:
    # migration parameters
    alpha12 = 0.0002
    alpha21 = 0.0002
if figure == 6:
    # migration parameters
    alpha12 = 0.0002
    alpha21 = 0.01
if figure == 7:
    # migration parameters
    alpha12 = 0.005
    alpha21 = 0.005 
if figure == 8:
    # migration parameters
    alpha12 = 0.1
    alpha21 = 0.1   
if figure == 9:
    # migration parameters
    alpha12 = 0.1
    alpha21 = 0.1   
    mu_s = 0.1
if figure == 10:
    # migration parameters
    alpha12 = 0.1
    alpha21 = 0.1   
    mu_s = 0.1
    T=30


# Time parameter
T = 50.0

# optimisation parameter
U_bar = 20000


# figure in tikz
tikz = False


# discretisation parameter
nt = int(300)


####################"
# resolution of the equilibrium
##############

m = GEKKO()
m.time = np.linspace(0,T,nt)

# Variables
P1 = m.Var(value=6000)
P2 = m.Var(value=6000)

# Vector to extract the final value
final_array = np.zeros(nt)
final_array[-1] = 1.0
final = m.Param(value=final_array)

# Equations
m.Equation(P1.dt() == b1*P1-mu11*P1-mu21*P1**2 + alpha21*P2 - alpha12*P1)
m.Equation(P2.dt() == b2*P2-mu12*P2-mu22*P2**2 + alpha12*P1 - alpha21*P2)

# Resolution
m.options.IMODE = 7
#m.options.NODES = 4
#m.options.MV_TYPE = 1
m.options.SOLVER = 2
#m.solver_options = ['max_iter 6000']
#m.options.OTOL = 0.1
print('begining of computations')
m.solve()
print('end of computations')

# Print results
P1_array=P1.value
P2_array=P2.value
t_array = np.arange(nt)*T/nt
P1T = float(P1_array[-1])
P2T = float(P2_array[-1])
print('P1(T): ' + str(P1T))
print('P2(T): ' + str(P2T))

# target
eps = P1T*0.01


#######################
# resolution of the optimal control
#####################

m = GEKKO()

# initial condition of the populations
P1_star =  float(P1_array[-1])
P2_star =  float(P2_array[-1])

# Time interval
m.time = np.linspace(0,T,nt)

# Control
u1 = m.MV(lb=0,ub=U_bar)
u1.STATUS = 1
u1.DCOST = 0.0
if test_case ==1:
	u2_array = np.zeros(nt)
	u2 = m.Param(value=u2_array)
else:
	u2 = m.MV(lb=0,ub=U_bar,value=np.zeros(nt))
	u2.STATUS = 1
	u2.DCOST = 0.0
print('bar U',U_bar)

# Variables
P1 = m.Var(value=P1T,lb=0.0)
P2 = m.Var(value=P2T,lb=0.0)
Ms1 = m.Var(value=0)
Ms2 = m.Var(value=0)
U = m.Var(value=0.0)#,lb=0,ub=C)

# Vector to extract the final value
final_array = np.zeros(nt)
final_array[-1] = 1.0
final = m.Param(value=final_array)
one_array = np.ones(nt)
one = m.Param(value=one_array)

# Equations
if test_case ==1:
	m.Equation(P1.dt() == b1*P1**2/(P1+gamma*Ms1)-mu11*P1-mu21*P1**2)
	m.Equation(Ms1.dt() == u1 - mu_s*Ms1)
	m.Equation(P2.dt() == b2*P2**2/(P2+gamma*Ms2)-mu12*P2-mu22*P2**2)
	m.Equation(Ms2.dt() == - mu_s*Ms2)
	m.Equation(U.dt() == u1 )
else:
	m.Equation(P1.dt() == b1*P1**2/(P1+gamma*Ms1)-mu11*P1-mu21*P1**2 + alpha21*P2 - alpha12*P1)
	m.Equation(Ms1.dt() == u1 - mu_s*Ms1+ omega*(alpha21*Ms2 - alpha12*Ms1))
	m.Equation(P2.dt() == b2*P2**2/(P2+gamma*Ms2)-mu12*P2-mu22*P2**2 + alpha12*P1 - alpha21*P2)
	m.Equation(Ms2.dt() == u2 - mu_s*Ms2+ omega*(alpha12*Ms1 - alpha21*Ms2))
	m.Equation(U.dt() == u1 + u2)
m.Equation(P1*final <= eps)

# Objective Function
m.Obj(U*final)

# Resolution
m.options.IMODE = 6
m.options.NODES = 3
#m.options.MV_TYPE = 1
m.options.SOLVER = 3
m.solver_options = ['max_iter 6000']
#m.options.OTOL = 0.1
print('begining of computations')
m.solve()
print('end of computations')

# Print results
P1_array=P1.value
Ms1_array=Ms1.value
P2_array=P2.value
Ms2_array=Ms2.value
u1_array=u1.value
u2_array=u2.value
U_array=U.value
t_array = np.arange(nt)*T/nt
J = float(P1_array[-1])
print('Objective = P1(T): ' + str(J))
print('int u1: ' + str(sum(u1_array)*float(T)/float(nt)))
print('int u2: ' + str(sum(u2_array)*float(T)/float(nt)))
print('int u: ' + str(sum(u1_array+u2_array)*float(T)/float(nt)))
print('int u: ' + str(U_array[-1]))


# Save the values for tikz
NN=5 #one value over 5
def write_to_file(f,A,B):
    f.write('\n')
    for i in range(len(A)):
        f.write('(')
        f.write(str(A[i]))
        f.write(',')
        f.write(str(B[i]))
        f.write(')')
        f.write('\n')
    f.write('\n')
    f.write('\n')
    return;
if tikz == True:
    f = open('simu_control_patch_test_case_{case}.txt'.format(case=str(test_case)),'w')
    f.write('Values of P1')
    f.write('\n')
    write_to_file(f,t_array[::NN],P1_array[::NN])
    f.write('\n')
    f.write('Values of P2')
    f.write('\n')
    write_to_file(f,t_array[::NN],P2_array[::NN])
    f.write('\n')
    f.write('Values of u1')
    f.write('\n')
    write_to_file(f,t_array[::NN],u1_array[::NN])
    f.write('\n')
    f.write('Values of u2')
    f.write('\n')
    write_to_file(f,t_array[::NN],u2_array[::NN])
    f.write('\n')
    f.close()

# Plot
plt.figure(1)
plt.plot(t_array,P1_array,'g-',linewidth=2,label=r'$P1$')
plt.plot(t_array,P2_array,'r-',linewidth=2,label=r'$P2$')
plt.legend(loc='best')
plt.xlabel('Time')
#plt.axis([0, T, min(E_array)-0.05*(max(E_array)-min(E_array)), max(E_array)+0.05*(max(E_array)-min(E_array))])
plt.savefig("simu_patch_control_P_T_{name0}_Ubar_{name1}_eps_{name2}_alpha12_{name3}_alpha21_{name4}.png".format(name0 = str(int(T)),name1 = str(U_bar),name2 = str(eps),name3 = str(alpha12),name4 = str(alpha21)))
plt.figure(2)
plt.plot(t_array,u1_array,'g-',linewidth=2,label=r'$u1$')
plt.plot(t_array,u2_array,'r-',linewidth=2,label=r'$u2$')
plt.legend(loc='best')
plt.xlabel('Time')
#plt.axis([0, T, min(F_array)-0.05*(max(F_array)-min(F_array)), max(F_array)+0.05*(max(F_array)-min(F_array))])
plt.savefig("simu_patch_control_u_T_{name0}_Ubar_{name1}_eps_{name2}_alpha12_{name3}_alpha21_{name4}.png".format(name0 = str(int(T)),name1 = str(U_bar),name2 = str(eps),name3 = str(alpha12),name4 = str(alpha21)))
plt.show()





