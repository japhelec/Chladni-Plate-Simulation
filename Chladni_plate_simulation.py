#-------------------------------1-------------------------------------
#import and parameter definition
#note that since we are simulating rope oscillation, we should set the rope elastic constant to be larger, in order to reproduce
#correct simulation
# A = oscillatopn amplitude, NN = number of particles, omega = oscillation angular frequency, N = square root of NN
# size = particle size, m = particle mass, k = bonding strength, d = original particle distance

from vpython import *
import numpy as np
A, NN, f = 0.4, 625, 5.725
size, m, k, d = 0.3, 0.1, 7500.0, 0.4
b = 2 #damping coefficient
N = int(sqrt(NN))
middle = int((N+1)/2)
array_m = middle-1
t, dt = 0, 0.001


#----------------------- 2 initialization ---------------------------
#scene and balls and springs settlement and initialiation
#set [m][n] as m+1 row n+1 column
scene = canvas(title='Chladni Pattern', width=800, height=800, background=vec(0.5,0.5,0), center = vec((N-1)*d/2, 0, 0))


# ----2-1. establish balls and springs
balls = []
for i in range(N):
    ball = [sphere(radius=size, color=color.red, pos=vector(2*j*d, 0, -2*i*d), v=vector(0,0,0)) for j in range(N)]
    balls.append(ball)
rsprings = []
csprings = []
'''
for i in range(N):
    springs = [helix(radius = size/2.0, thickness = d/15.0, pos=vector(j*d, 0, -i*d), axis=vector(d,0,0)) for j in range(N-1)]
    rsprings.append(springs)

for i in range(N-1):
    springs = [helix(radius = size/2.0, thickness = d/15.0, pos=vector(j*d, 0, -i*d), axis=vector(0,0,-d)) for j in range(N)]
    csprings.append(springs)
'''


# ---2-2. ballsarray pos velocity 3D
ballspos = []
ballsv = []
for i in range(N):
    aux_pos = []
    aux_v = []
    for j in range(N):
        aux_pos.append([2*j*d, 0, -2*i*d])
        aux_v.append([0,0,0])
    ballspos.append(aux_pos)
    ballsv.append(aux_v)

ballspos = np.array(ballspos, dtype = float)
ballsv = np.array(ballsv, dtype = float)  #transform list to array



ballsamp = []
ballsramp = []
for i in range(N):
    aux_amp = []
    aux_ramp = []
    for j in range(N):
        aux_amp.append(0)
        aux_ramp.append(0)
    ballsamp.append(aux_amp)
    ballsramp.append(aux_ramp)
ballsamp = np.array(ballsamp, dtype = float)
ballsramp = np.array(ballsramp, dtype = float)


# ---2-3. springs :
#important parameters:
#  rspringslen_o:  row springs length at beginning 1D
#  cspringslen_o:
#  rspringsaxis :  row springs axis 3D
#  cspringsaxis :
#  cspringslen  :  row springs length now 1D
#  cspringslen  :
#  rspringsacc:  row springs acceleration now 3D
#  cspringsacc: 


# row springs original length and force
rspringslen_o = [] 
cspringslen_o = [] 
rspringsacc = []
cspringsacc = []

for i in range(N):
    raux_leno = []
    raux_acc = []
    for j in range(N-1):
        raux_leno.append(d)
        raux_acc.append([0,0,0])
    rspringslen_o.append(raux_leno)
    rspringsacc.append(raux_acc)

for i in range(N-1):
    caux_leno = []
    caux_acc = []
    for j in range(N):
        caux_leno.append(d)
        caux_acc.append([0,0,0])
    cspringslen_o.append(caux_leno)
    cspringsacc.append(caux_acc)

rspringslen_o = np.array(rspringslen_o, dtype = float)
cspringslen_o = np.array(cspringslen_o, dtype = float)
rspringsacc = np.array(rspringsacc, dtype = float)
cspringsacc = np.array(cspringsacc, dtype = float)

#spring axis and acceleration
rspringsaxis = ballspos[:,1:,:]-ballspos[:,:-1,:]
raux1 = (rspringsaxis)**2
rspringslen = np.sqrt(raux1[:,:,0]+raux1[:,:,1]+raux1[:,:,2])
raux2 = k*(rspringslen-rspringslen_o)/rspringslen/m*dt
rspringsacc[:,:,0] = raux2*rspringsaxis[:,:,0]

rspringsacc[:,:,1] = raux2*rspringsaxis[:,:,1]
rspringsacc[:,:,2] = raux2*rspringsaxis[:,:,2]


cspringsaxis = ballspos[1:,:,:]-ballspos[:-1,:,:]
caux1 = (cspringsaxis)**2
cspringslen = np.sqrt(caux1[:,:,0]+caux1[:,:,1]+caux1[:,:,2])
caux2 = k*(cspringslen-cspringslen_o)/cspringslen/m*dt
cspringsacc[:,:,0] = caux2*cspringsaxis[:,:,0]
cspringsacc[:,:,1] = caux2*cspringsaxis[:,:,1]
cspringsacc[:,:,2] = caux2*cspringsaxis[:,:,2]
#print(cspringsacc)

n2 = 0
n = 50
#n2 = 10
#-----------------------------------3 run the simulation---------------------------
op = open('2Dfreq.txt', 'w')
while t < 20:
    print(t)
    t += dt

    # let middle ball oscillates sinusoidally
    ballspos[array_m][array_m][1] = A * sin(2*pi*f * t)

    #refresh springs:
    rspringsaxis = ballspos[:,1:,:]-ballspos[:,:-1,:]
    raux1 = (rspringsaxis)**2
    rspringslen = np.sqrt(raux1[:,:,0]+raux1[:,:,1]+raux1[:,:,2])
    raux2 = k*(rspringslen-rspringslen_o)/rspringslen/m*dt
    
    rspringsacc[:,:,0] = raux2*rspringsaxis[:,:,0]
    rspringsacc[:,:,1] = raux2*rspringsaxis[:,:,1]
    rspringsacc[:,:,2] = raux2*rspringsaxis[:,:,2]


    cspringsaxis = ballspos[1:,:,:]-ballspos[:-1,:,:]
    caux1 = (cspringsaxis)**2
    cspringslen = np.sqrt(caux1[:,:,0]+caux1[:,:,1]+caux1[:,:,2])
    caux2 = k*(cspringslen-cspringslen_o)/cspringslen/m*dt
    cspringsacc[:,:,0] = caux2*cspringsaxis[:,:,0]
    cspringsacc[:,:,1] = caux2*cspringsaxis[:,:,1]
    cspringsacc[:,:,2] = caux2*cspringsaxis[:,:,2]
    


    #refresh balls:
    
    ballsv[1:-1,1:-1,:] += rspringsacc[1:-1,1:,:] - rspringsacc[1:-1,0:-1,:] + cspringsacc[1:,1:-1,:] - cspringsacc[0:-1,1:-1,:]   -b*ballsv[1:-1,1:-1,:]*dt 
    ballsv[array_m,array_m,:] = [0,0,0]    
    ballspos[:,:,:] += ballsv[:,:,:]*dt
    ballsamp = np.fmax(ballsamp,ballspos[:,:,1])
    ballsramp[1:-1,1:-1] = ballsamp[1:-1,1:-1]/np.max(ballsamp[1:-1,1:-1])

    op.write(str(ballspos[array_m,array_m,1]) + '\n')

        
    # #present balls:
    # if t > n*dt:
    #     for i in range(N):
    #         for j in range(N):
    #             balls[i][j].pos = vec(ballspos[i,j,0], 0, ballspos[i,j,2])
    #             balls[i][j].color = color.gray(ballsramp[i,j])
    #     n += 50

op.close()
quit()