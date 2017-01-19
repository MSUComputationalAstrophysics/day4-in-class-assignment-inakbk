#in class 4

from pylab import *

#----------------------------------------------------------------------------------
#analytical solutions
def x(t, x_m, omega):
	return x_m*sin(omega*t)
def v(t, x_m, omega):
	return x_m*cos(omega*t)

#calculating total energy:
def E_tot(x, v, k, m):
	return 0.5*k*x*x + 0.5*m*v*v
#calculating relative change in energy:
def epsilon(E_start, E_end):
	return abs(E_end - E_start)/E_start

#----------------------------------------------------------------------------------
#ODE methods:

#the acceleration (in the pre-class assignment I used f in stead, but a = f/m):
def a(x, k, m): 
	return -k*x/m

def Euler(x_0, v_0, a, len_t, dt, m, k):
	x_num = zeros(len_t)
	v_num = zeros(len_t)
	x_num[0] = x_0
	v_num[0] = v_0
	for i in range(len_t-1):
		x_num[i+1] = x_num[i] + dt*v_num[i]
		v_num[i+1] = v_num[i] + dt*a(x_num[i],k, m)
	#calculating start and end energy
	E_start =  E_tot(x_num[0], v_num[0], k, m)
	E_end =  E_tot(x_num[-1], v_num[-1], k, m)
	return x_num, E_start, E_end

def RK4(x_0, v_0, a, len_t, dt, m, k):
	x_num = zeros(len_t)
	v_num = zeros(len_t)
	x_num[0] = x_0
	v_num[0] = v_0
	for i in range(len_t-1):
		x_1 = x_num[i]
		v_1 = v_num[i]
		a_1 = a(x_1, k, m) #a flops
		#-------------------
		x_2 = x_1 + v_1*dt/2 #4flops
		v_2 = v_1 + a_1*dt/2 #4flops
		a_2 = a(x_2, k, m) #a flops
		#-------------------
		x_3 = x_1 + v_2*dt/2 #4flops
		v_3 = v_1 + a_2*dt/2 #4flops
		a_3 = a(x_3, k, m) #a flops
		#-------------------
		x_4 = x_1 + v_3*dt #3 flops
		v_4 = v_1 + a_3*dt #3 flops
		a_4 = a(x_4, k, m) #a flops
		#-------------------
		x_num[i+1] = x_1 + (1/6.)*(v_1 + 2*v_2 + 2*v_3 + v_4)*dt #10 flops
		v_num[i+1] = v_1 + (1/6.)*(a_1 + 2*a_2 + 2*a_3 + a_4)*dt #10 flops
	#calculating start and end energy
	E_start =  E_tot(x_num[0], v_num[0], k, m)
	E_end =  E_tot(x_num[-1], v_num[-1], k, m)
	return x_num, E_start, E_end

#----------------------------------------------------------------------------------
#constants:
m = 1 #mass
k = 1 #spring constant
omega = 1 #sqrt(k/m) =1 in this case
x_m = 1

#initial conditions and end time:
t_0 = 0 
x_0 = 0 
v_0 = 1
t_end = 4 #scaling the time with pi

#----------------------------------------------------------------------------------
#doing all the timesteps in one go
dt_list = linspace(0.7, 0.0001,20) #[0.1, 0.01, 0.001] 
col_list = ['b', 'g', 'k']

#storing values of the energy calculated after the integration
E_start_euler = zeros(len(dt_list))
E_end_euler = zeros(len(dt_list))
E_start_RK4 = zeros(len(dt_list))
E_end_RK4 = zeros(len(dt_list))

#---------------------------------------------------------------------------------
#calculating for all timesteps for euler and RK4:
#figure(1)
for j in range(len(dt_list)): 
	N = (t_end - t_0)/dt_list[j]
	t = linspace(t_0, t_end, N+1)
	#calculating the position and start/end energy:
	x_RK4, E_start_RK4[j], E_end_RK4[j] = RK4(x_0, v_0, a, len(t), dt_list[j], m, k)
	x_euler, E_start_euler[j], E_end_euler[j] = Euler(x_0, v_0, a, len(t), dt_list[j], m, k)
	#plot(t, x_RK4, col_list[j])
	#plot(t, x_euler, col_list[j])
#analytical:
#plot(t, x(t, x_m, omega), 'r--') 
#legend(['dt=0.1', 'dt=0.01', 'dt=0.001', 'analytical'])
#title('Position versus time with the Euler method')
#xlabel('t/pi')
#ylabel('x')

print 'Difference in relative change in energy between Euler and RK4:' 
print abs(epsilon(E_start_euler, E_end_euler) - epsilon(E_start_RK4, E_end_RK4))
print 'for the timesteps [0.1, 0.01, 0.001]'

#plotting the change in total error:
figure(2)
loglog(dt_list, epsilon(E_start_euler, E_end_euler) , 'bo-') 
loglog(dt_list, epsilon(E_start_RK4, E_end_RK4) , 'ro-') 
loglog(dt_list, 1e-2*ones(len(dt_list)), 'k')

legend(['Euler', 'RK4', 'Accuracy wanted'], loc='center left')
title('Relative change in total energy versus timestep')
xlabel('dt/pi')
ylabel('epsilon')

show()
















