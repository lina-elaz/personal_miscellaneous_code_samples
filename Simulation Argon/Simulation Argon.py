#*****************************************************************************80
#
#  MD_3D is the main program for simulating molecular dynamcis.
#
#  Discussion:
#    Basic (microcanonical) MD program for N point particles that
#    interact via a (shifted & truncated) purely repulsive Lennard-
#    Jones potential.

#
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    08 December 2015
#
#  Author:
#
#    MEHDI AYOUZ LGPM CENTRALESUPELEC
#    mehdi.ayouz@centralesupelec.fr
#    Bat DUMAS C421. 0141131603
#
#  Parameters:
#   sigma, epsilon, m, dt, N, maxstep, step, Nrdf, rmax (see lecture notes)
#
#
#
from __future__ import division
from math import *
import numpy as np
from numpy import *
import cmath as c
from cmath import *
from pylab import *
import pylab as pl
from numpy import linalg as LA
from time import clock
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import random as rd

#*****************************************************************************#
class Get_V_or_F(object):

	def __init__(self,rij):
		self.rij = rij

	def V_Lennard_Jones(self):

		XRij = (4*epsilon*( (sigma/self.rij)**12 - (sigma/self.rij)**6)+epsilon)

		return XRij

	def F_Lennard_Jones(self):

		XRij = 48.0e0*epsilon*((sigma/self.rij)**12/self.rij - 0.5*(sigma/self.rij)**6/self.rij )

		return XRij

#*****************************************************************************#

def ReplaceInBox(X,Y,Z,L): #MY PART
    Xtemp = X + L/2
    Ytemp = Y + L/2
    Ztemp = Z + L/2
    Xtemp = np.mod(Xtemp, L)
    Ytemp = np.mod(Ytemp, L)
    Ztemp = np.mod(Ztemp, L)
    X = Xtemp - L/2
    Y = Ytemp - L/2
    Z = Ztemp - L/2
    return X, Y, Z

# this is the application of the PBC (periodic boundaries condition)

#*****************************************************************************#
def SpecialDistance(X1, X2, Y1, Y2, Z1, Z2, L): #MY PART
    if np.absolute(X2 - X1) > np.absolute(X2 + L - X1):
        R21X = X2 + L - X1
    elif np.absolute(X2 - X1) > np.absolute(X2 - L - X1):
        R21X = X2 - L - X1
    else:
        R21X = X2 - X1

    if np.absolute(Y2 - Y1) > np.absolute(Y2 + L - Y1):
        R21Y = Y2 + L - Y1
    elif np.absolute(Y2 - Y1) > np.absolute(Y2 - L - Y1):
        R21Y = Y2 - L - Y1
    else:
        R21Y = Y2 - Y1

    if np.absolute(Z2 - Z1) > np.absolute(Z2 + L - Z1):
        R21Z = Z2 + L - Z1
    elif np.absolute(Z2 - Z1) > np.absolute(Z2 - L - Z1):
        R21Z = Z2 - L - Z1
    else:
        R21Z = Z2 - Z1

    return R21X, R21Y, R21Z
    
    # PBC implies special distance, the minimal distance is taken into account 

def timestamp ( ):
#*****************************************************************************
#
## TIMESTAMP prints the date as a timestamp.
#
#  Licensing:
#
#   This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#   11 October 2015
#
#  Author:
#
#   Mehdi Ayouz
#
#  Parameters:
#
#   None
#
    import time

    t = time.time ( )
    print(time.ctime ( t ))

    return None

#   Generate the initial conformation (positions and velocities).
#   Generate initial positions on sites of a simple cubic lattice
#       (plus additional small random displacements):
#       NOTE: Coordinates in simulation box are -Lx/2 <= rX < +Lx/2, etc.
#             (Here even -Lx/2 < rX, etc., holds.)

def Get_Velocity_Distibution(V,Nv=101):
#*****************************************************************************#

    maxval=max(V)
    minval=min(V)

#   Plot velocity distribution
    figure(5)
    hist(V, Nv, [minval, maxval], normed=1, alpha=0.5, label="Velocity Distribution")
    legend()
    xlabel('V')
    ylabel('f(V)')
    draw()
    show()


    return None
#*****************************************************************************#
def Get_Initial_Conformation(N,Lx,rc,m,Etot,Lxbox,Tfac, sampling): #MY PART

    Nlattice= int(N**(1.0/3.0)) + 1
    a=Lx/Nlattice
    iplaced = 0
    
    # iplaced count the number of particles are in place

    X = np.zeros(N,float)
    Y = np.zeros(N,float)
    Z = np.zeros(N,float)

    Vec_particles=np.zeros(N,int)
    Lybox=Lxbox
    Lzbox=Lybox
    Ly=Lx
    Lz=Ly

    print('a= %g' % a)

#   a is the distance from -L/2  to +L/2 with small gaussian displacemet
    for i in range(0,Nlattice):
        for j in range(0,Nlattice):
            for k in range(0,Nlattice):
                if (iplaced < N):
                    X[iplaced]= -Lx/2 + float(i)*a + rd.gauss(0,0.005)
                    Y[iplaced]= -Lx/2 + float(j)*a + rd.gauss(0,0.005)
                    Z[iplaced]= -Lx/2 + float(k)*a + rd.gauss(0,0.005)

                    X[iplaced], Y[iplaced], Z[iplaced] = ReplaceInBox(X[iplaced], Y[iplaced], Z[iplaced], Lx)

                    iplaced += 1




#   Verify that ALL the particles of the system are in place
    print('iplaced=%g' % iplaced)

#    initialize random particle velocities .shifting
#    and re-scaling  done later:

    VX = np.zeros(N, float)
    VY = np.zeros(N, float)
    VZ = np.zeros(N, float)
    V  = np.zeros(N, float)
    for i in range(0,N):
        if (sampling == 'Uniform'):
            def RdGen():
                return rd.random() * 2 - 1
        elif (sampling == 'Gaussian'):
            def RdGen():
                return rd.gauss(0,1)

        else:
            print ('Sampling must be Gaussian or Uniform; Please check input variable in function InitVelocity')
        VX[i] = RdGen()
        VY[i] = RdGen()
        VZ[i] = RdGen()


        V[i]=sqrt(VX[i]**2+VY[i]**2+VZ[i]**2)

#    Correction of generated velocities in order to have:
#    center-of-mass velocity =0, and Epot+Ekin = to the set total energy value Etot


# Calculate current center of mass velocity
    VCoMX = sum(VX)
    VCoMY = sum(VY)
    VCoMZ = sum(VZ)

#  As  the center of mass velocity= 0  Correct all velocities by VcoMX, VcoMY and VcoMZ
    VX    = VX - VCoMX/float(N)  # Removing the center of mass velocity
    VY    = VY - VCoMY/float(N)
    VZ    = VZ - VCoMZ/float(N)




# Calculate potential energy per particle by summing over
#     ALL pairs of particles.

    Epot = 0.0
    for j in range(0,N):
        for i in range(0,N):
            if i != j:
                Rijx, Rijy, Rijz = SpecialDistance(X[i], X[j], Y[i], Y[j], Z[i], Z[j], Lx)
                Rij = np.sqrt(Rijx**2 + Rijy**2 + Rijz**2)
                if(Rij < rc):
                    Force = Get_V_or_F(Rij)
                    Epot += Force.V_Lennard_Jones()

    Epot=Epot/N  # this is potential energy per particle

    if(Etot<Epot):
        print('Attention : Epot=%g'%Epot)
        print('Etot = %g'%Etot)


#  In order to have Ekin = Etots - Epot we must Re-scale velocities by ScaleVelocity :


    EkinTemp =  0.5e0*m*sum(VX**2 + VY**2 + VZ**2)/float(N) # kinetic energy per particle, so are Epot and Etot 
    ScalVelocity = np.sqrt( (Etot-Epot)/EkinTemp)
    VX = ScalVelocity*VX
    VY = ScalVelocity*VY
    VZ = ScalVelocity*VZ


#    Recalculate kinetic energy from shifted and re-scaled velocities:

    V  = np.sqrt(VX**2 + VY**2 + VZ**2)
    Ekin =  0.5e0*m*sum(VX**2 + VY**2 + VZ**2)/float(N) 
    # kinetic energy per particle

    temp = 2*Ekin* N /(3*N - 3)          #  temperature thanks to Ekin

    Etot = Epot + Ekin                   #  per particle


#   Check out the velocity distibution
    Get_Velocity_Distibution(V,N/2 + 1)

#   Plot initial conformation (xi,yi,zi)
    figure(6) #TBD


    return X,Y,Z,VX,VY,VZ,Ekin,Epot,Etot,temp
#*****************************************************************************#
def Get_Forces(N,rc,Lx,X,Y,Z): #MY PART

#   FX , FY and FZ are tables full of zeros
    FX = np.zeros(N, float)
    FY = np.zeros(N, float)
    FZ = np.zeros(N, float)
    Fij = np.zeros((N, N), float)

    Ly=Lx;Lz=Ly
    W=0.
    Epot=0.

    for j in range(0,N):
        for i in range(0,N):
            if(i != j):
                Rijx, Rijy, Rijz = SpecialDistance(X[i], X[j], Y[i], Y[j], Z[i], Z[j], Lx)
                Rij = np.sqrt(Rijx**2 + Rijy**2 + Rijz**2)
                if(Rij < rc):
                    Force = Get_V_or_F(Rij)
                    Fij[i,j] = Force.F_Lennard_Jones() # square matrix of dimension N and each term is proportionnal to the force from j to i
                    W += Fij[i,j]
                    Fij = Fij/Rij**2
                    FX[i] += -Fij[i,j]*Rijx
                    FY[i] += -Fij[i,j]*Rijy
                    FZ[i] += -Fij[i,j]*Rijz

    return FX,FY,FZ,W 
    # we have the forces applied on each particles according to Ox, Oy and Oz

#*****************************************************************************#
def Get_Velocity_Integration(N,VX,VY,VZ,FX,FY,FZ,dtH): #dtH= dt/2    #MY PART

# Use of Verlet Scheme
    V = np.zeros(N, float)
    for i in range(0, N):
        VX[i] = VX[i] + FX[i]*dtH/m
        VY[i] = VY[i] + FY[i]*dtH/m
        VZ[i] = VZ[i] + FZ[i]*dtH/m
        V[i] = np.sqrt(VX[i]**2 + VY[i]**2 + VZ[i]**2)

    return VX,VY,VZ,V # we have the velocity of each particle
#*****************************************************************************#
def Get_Position_Integration(N,X,Y,Z,VX,VY,VZ,dt): #MY PART
    
    # Use of verlet scheme
    for i in range(0, N):
        X[i] = X[i] + VX[i]*dt
        Y[i] = Y[i] + VY[i]*dt
        Z[i] = Z[i] + VZ[i]*dt

    return X,Y,Z
#*****************************************************************************#
def Get_Rdial_Distribution_Function_Sample(N,Lx,rmax,Bw,X,Y,Z):

    Ly=Lx;Lz=Ly

    grdf = np.zeros(Nrdf, float)

    for ip in range(N-1):
        Xip=X[ip]
        Yip=Y[ip]
        Zip=Z[ip]

        for jp in range(ip+1,N):
            Xij=Xip-X[jp]
            Yij=Yip-Y[jp]
            Zij=Zip-Z[jp]

            Xij=Xij - Lx*round(Xij/Lx)
            Yij=Yij - Ly*round(Yij/Ly)
            Zij=Zij - Lz*round(Zij/Lz)

            Rij=sqrt(Xij**2+Yij**2+Zij**2)

            if(Rij<=rmax):
                irdf=int(Rij/Bw)
                grdf[irdf-1]=grdf[irdf-1]+2

    return grdf

def Get_Rdial_Distribution_Function(Rdffac,Bw,Nrdf,grdf):

    Rrdf=np.zeros(Nrdf,float)
    g=np.zeros(Nrdf,float)
    d=np.zeros(Nrdf, float)

    for i in range(1,Nrdf):
        Rrdf[i]=(float(i+1)-0.5)*Bw
        Vshell=((i+1)*Bw)**3-((i)*Bw)**3
        g[i]=grdf[i]*Rdffac/Vshell
        d[i]=Bw*i

#   Plot g(r)
    figure(2)
    plot(d,g,'-b', label="Radial Distribution")


    return None
#*****************************************************************************#
#  Calculate instantaneous values of
#  - modulus of center-of-mass velocity   Vcm
#  - kinetic energy per bead             Ekin
#  - potential energy per bead           Epot
#  - total energy per bead               Etot
#  - (kinetic) temperature               temp
#  - pressure                           press
#*****************************************************************************#
def Get_Observables(N,rc,Tfac,Pfac,X,Y,Z,VX,VY,VZ, W):  #MY PART
    VXcm= 0.0
    VYcm= 0.0
    VZcm= 0.0
    Ekin = 0.0
    Vcm=0.

    VXcm = sum(VX)/float(N)
    VYcm = sum(VY)/float(N)
    VZcm = sum(VZ)/float(N)

    Vcm = np.sqrt(VXcm**2 + VYcm**2 + VZcm**2) # modulus of Vcm vector

    Ekin =  0.5e0*m*sum(VX**2 + VY**2 + VZ**2)/float(N) # kinetic energy per particle
    temp = 2*Ekin* N /(3*N - 3)  #  temperature thanks to Ekin
    press = rho*temp+Pfac*W # pressure

    Epot = 0.0
    for i in range(0,N):
        for j in range(0,N):
            if(i != j):
                Rijx = X[i] - X[j]
                Rijy = Y[i] - Y[j]
                Rijz = Z[i] - Z[j]
                Rij = np.sqrt(Rijx**2 + Rijy**2 + Rijz**2)
                if(Rij<rc):
                    Epot += Get_V_or_F(Rij).V_Lennard_Jones() 
#  Lennard Jones interaction as described in Mr Ayouz' Report

    Epot=Epot/N  # potential energy per particle



    Etot = Epot + Ekin              # total energy per particle
    return Ekin,Epot,Etot,temp,press,Vcm
#*****************************************************************************#
def Get_Self_Diffusivity(N,ndata,ntime,step,dt,VecR_msd):

    print("read all the data for mean sqaure displacments calculations ")
#   Number of origins is half number of time steps
    norigin = int(ntime/2)

#   Minimum number of intervals to contribute to diffusivity
    nmin = 50

#   Maximum number of intervals to contribute to diffusivity
    nmax = norigin

    if(nmin>nmax):
        print("We have a problem in computing Self-Diffusivity ")
        print("nmin = %g" % nmin, " nmax = %g" % nmax)

#   Store mean square displacements in xmsd
    time_vec=np.zeros(norigin,float)
    Rmsd=np.zeros((norigin,3),float)
    slope=np.zeros(3,float)
    slopesd=np.zeros(3,float)
    yinter=np.zeros(3,float)
    yintersd=np.zeros(3,float)
    Dav=np.zeros(4,float)
    Dsd=np.zeros(4,float)

#    for i in range(ndata):
#        print(VecR_msd[i,0],VecR_msd[i,1],VecR_msd[i,2])


    for i in range(norigin):
        time_vec[i] = ((i+1)*step)*dt

    for i in range(N):
        for j in range(norigin):
            jstart=j*N+i

            for k in range(nmin,nmax):
                kend=jstart+k*N

                Rmsd[k,0]=Rmsd[k,0]+(VecR_msd[kend,0]-VecR_msd[jstart,0])**2
                Rmsd[k,1]=Rmsd[k,1]+(VecR_msd[kend,1]-VecR_msd[jstart,1])**2
                Rmsd[k,2]=Rmsd[k,2]+(VecR_msd[kend,2]-VecR_msd[jstart,2])**2

    Rmsd=Rmsd/(N*norigin)

#   Perform a linear least squares regression
    tmp_Rmsd=np.zeros((nmax-nmin),float)
    tmp_time_vec=np.zeros((nmax-nmin),float)
    for i in range(3):
        for j in range(nmin,nmax):
            tmp_Rmsd[j-nmin]=Rmsd[j,i]
            tmp_time_vec[j-nmin]=time_vec[j]

        tmp_slope,tmp_slopesd,tmp_yinter,tmp_yintersd=Get_Linear_Least_Square_Regression(nmax-nmin,tmp_time_vec,tmp_Rmsd)
        slope[i]=tmp_slope
        slopesd[i]=tmp_slopesd
        yinter[i]=tmp_yinter
        yintersd[i]=tmp_yintersd

    charname=['x','y','z','avg']

    for i in range(3):
        print("%s "%charname[i],"Slope=%g" %slope[i], "y-intercept=%g"%yinter[i])
        Dav[i]= 0.50*slope[i]*msdfac  # convert to m^2/sec
        Dsd[i]= 0.50*slopesd[i]*msdfac # convert to m^2/sec

    Dav[3] = sum(Dav[0:2])/3.0

#   Standard deviation of average diffusivity
    term1 = 3.0*(Dav[0]**2 + Dav[1]**2 + Dav[2]**2)
    Dsd[3] = sqrt( (term1 - Dav[3]*Dav[3]*9.0) /6.0 )

    for i in range(4):
        print("%s "%charname[i],"Diffusivity avg=%g" %Dav[i], "standar deviation =%g" %Dsd[i],"in m^2/sec")

#   Write xmsd vs time data for plotting
    figure(4)
    plot(time_vec, Rmsd[:,0],'-b',time_vec, Rmsd[:,1],'-r',time_vec, Rmsd[:,2],'-g',label="MSD")
    title("MSD vs time")
    legend()
    xlabel("time")
    ylabel("MSD")


    return None
#*****************************************************************************#
def Get_Linear_Least_Square_Regression(n, x, y):

    xn = float(n)
    xavg = sum(x)/xn
    yavg = sum(y)/xn
    sumxy = 0.0
    sumxx = 0.0
    sumx2 = 0.0

    for i in range(n):
        sumxy = sumxy + (x[i] - xavg)*(y[i] - yavg)
        sumxx = sumxx + (x[i] - xavg)*(x[i] - xavg)
        sumx2 = sumx2 + x[i]*x[i]

    tmp_slope = sumxy/sumxx
    tmp_yinter = yavg - tmp_slope*xavg
    sse = 0.0

    for i in range(n):
        sse = sse + (y[i] - tmp_slope*x[i] -tmp_yinter)**2.0

    sig2 = sse/float(n-2)
    tmp_slopesd = sqrt(sig2/sumxx)
    tmp_yintersd = sqrt(sig2/float(n)*sumx2/sumxx)

    return  tmp_slope,tmp_slopesd,tmp_yinter,tmp_yintersd
#*****************************************************************************#
if __name__ == "__main__":

    timestamp ( )
    print('')
    print ('Molecular Dynamics simualtion')
    print ('Repulsive Lennard-Jones potential.')

    wtime1 = clock ( )

#   Open datafile
    f0 = open('Params.in', 'r')
#   Read datafile
    for line in f0:
        DumList = line.split()
        tag     = DumList[0]
        if (tag == 'N'):
            N       = int(DumList[1]) # total bead number
        elif (tag == 'rho'):
            rho     = float(DumList[1]) # global bead number density
        elif (tag == 'Etot'):
            Etot    = float(DumList[1])  # set total energy per particle.
        elif (tag == 'm'):
            m       = float(DumList[1]) #  mass
        elif (tag == 'dt'):
            dt      = float(DumList[1]) # integration time step
        elif (tag == 'maxstep'):
            maxstep = int(DumList[1]) # total no. of time steps (-> tmax = maxstep*dt)
        elif (tag == 'step'):
            step    = int(DumList[1]) # write output after each block of incstp time steps
        elif (tag == 'rmax'):
            rmax    = float(DumList[1]) # maximum radius for RDF sampling
        elif (tag == 'Nrdf'):
            Nrdf    = int(DumList[1]) # no. of bins for RDF sampling
        elif (tag == 'sigma'):
            sigma   = float(DumList[1]) # Lennard-Jones potential lenght
        elif (tag == 'epsilon'):
            epsilon = float(DumList[1]) # Lennard-Jones potential depth
#   Declare variables


    dtH=dt/2

    if(rho<0.01):
        print('Attention : rho too small %g' % rho)
    sampling = "Gaussian"
    blend = False
    Volume=N/rho
    L=Volume**(1./3.)
    Lx=L
    Ly=Lx
    Lz=Ly
    Lbox=L/2.
    Lxbox=Lbox
    Lybox=Lbox    # Cooridnate within the simulation box -Lx/2<=x<=Lx/2
    Lzbox=Lbox
    rc=2**(1./6.)*sigma




    timefac=2.17e-12    # in second
    spacefac=3.405e-10  # in m
    massfac=6.69e-26    # in kg=39.948 uma
    energyfac=0.01032   # in eV
    
    msdfac=spacefac**2/timefac

#   Number of times represented in data for data mean square displacments
    ntime = int(maxstep/step)+1

#   Number of rows of data mean square displacments
    ndata = N*ntime
    VecR_msd=np.zeros((ndata,3),float)

    if(maxstep % step != 0 ):
         print('MD_3D: incstp<>maxstp')
    tmax=dt*float(maxstep)
    Rfac= sqrt(3.0)
    Pfac= 1.0/(3.0*Volume)
    Tfac= 1.0/(3.0*N - 3.0)

#   Constants for RDF calculation:
    Rdffac = 0.75/(rho*N*pi)
    Bw = rmax/float(Nrdf-1)

    print('Volume=%g'% Volume)
    print('L=%g'% L)
    print('Lxbox=Lybox=Lzbox=%g'% Lxbox)
    print('rc=%g'% rc)

    nconf = 0;imsd=0

#   Generate initial particle conformation:
    X,Y,Z,VX,VY,VZ,Ekin,Epot,Etot,temp = Get_Initial_Conformation(N,Lx,rc,m,Etot,Lxbox,Tfac, sampling)


#   Compute forces from initial particle conformation:
    FX,FY,FZ,W = Get_Forces(N,rc,Lx,X,Y,Z)

    press=rho*temp+Pfac*W


    print ('    Epot = %g'  % Epot)
    print ('    Ekin = %g'  % Ekin)
    print ('    Etot = %g'  % Etot)
    print ('    temp = %g'  % temp)
    print ('    press = %g'  % press)

#   Zero averages:

    tempav  = 0.0     #  average (kinetic) temperature
    pressav = 0.0     #  average pressure
    Epotav  = 0.0     #  average potential energy
    Ekinav  = 0.0     #  average kinetic energy
    Etotav  = 0.0     #  average total energy

    Time = np.zeros(ntime - 1, float)
    Temp = np.zeros(ntime - 1, float)
    P = np.zeros(ntime - 1, float)
    K = np.zeros(ntime - 1, float)
    Vpot = np.zeros(ntime - 1, float)
    E = np.zeros(ntime - 1, float)

#   Loop over integration steps:
    for it in range(0, maxstep):

        print ('it= %d' % it)

#       Propagate velocities over dt/2:     # velocity-Verlet integratio starts here
        VX,VY,VZ,V = Get_Velocity_Integration(N,VX,VY,VZ,FX,FY,FZ,dtH)
#       Propagate positions over dt:
        X, Y, Z = Get_Position_Integration(N,X,Y,Z,VX,VY,VZ,dt)
        for k in range(0, N):
            X[k], Y[k], Z[k] = ReplaceInBox(X[k], Y[k], Z[k], L)

#       Compute forces:
        FX, FY, FZ, W = Get_Forces(N,rc,Lx,X,Y,Z)

#       Propagate velocities over dt/2:
        VX,VY,VZ,V = Get_Velocity_Integration(N,VX,VY,VZ,FX,FY,FZ,dtH)

#       Every "incstep" steps, sample data, and produce output:
        if(it%step == 0):

#           Increment counter for conformations used for data sampling:
            nconf = nconf + 1

#           Compute center-of-mass velocity, temperature, pressure,
#           and energies (kinetic, potential, total) per bead:
            Ekin,Epot,Etot,temp,press,Vcm = Get_Observables(N,rc,Tfac,Pfac,X,Y,Z,VX,VY,VZ, W)

#           Sample radial distribution function (rdf).
            grdf=Get_Rdial_Distribution_Function_Sample(N,Lx,rmax,Bw,X,Y,Z)


#           Normalize averages and write them to datfile:
            tempav  += temp
            pressav += press
            Epotav  += Epot
            Ekinav  += Ekin
            Etotav  += Etot

#           Save Observables time, Temp, P, K,  Vpot and E, to be plotted
            Time[it]=it*dt
            Temp[it]=temp
            P[it]=press
            K[it]=Ekin
            Vpot[it]=Epot
            E[it]=Etot

#           Store all R data for mean square displacements
            for ip in range(N):
                VecR_msd[ip+it*N/step,0]=X[ip] #ip*ntime+it ip+it
                VecR_msd[ip+it*N/step,1]=Y[ip]
                VecR_msd[ip+it*N/step,2]=Z[ip]
                imsd=imsd+1

        print('imsd=%d' %(imsd))
        print('Np*ntime=%d' %(N*ntime))

#   Compute self-diffusivity
    Get_Self_Diffusivity(N,ndata,ntime,step,dt,VecR_msd)

    print ('    Epotav = %g'  % (Epotav/nconf))
    print ('    Ekinav = %g'  % (Ekinav/nconf))
    print ('    Etotav = %g'  % (Etotav/nconf))
    print ('    tempav = %g'  % (tempav/nconf))
    print ('    pressav = %g'  % (pressav/nconf))
    print ('    Vcm = %g'%Vcm)

#   Plot Observables versus time
    figure(1)
    plot(Time, Temp, label="Temperature")
    title("Temperature vs time")
    legend()

    figure(10)
    plot(Time, P, label="Pressure")
    title("Pressure vs time")
    legend()


    figure(11)
    plot(Time, K, label="Kinetic Energy")
    title("EKin vs time")
    legend()


    figure(12)
    plot(Time, Vpot, label="Potential Energy")
    title("Epot vs time")
    legend()


    figure(13)
    plot(Time, E, label="Total Energy" )
    title("Etot vs time")
    legend()



#   Compute RDF and write data to file:
    Rdffac=Rdffac/nconf
    Get_Rdial_Distribution_Function(Rdffac,Bw,Nrdf,grdf)

#   Get velocity distibution function
    Get_Velocity_Distibution(V,Nv=101)

#   Plot particle positions (xi,yi,zi)
    figure(7)
    ax=plt.axes(projection='3d')
    ax.plot(X,Y,Z,'o',label="Particles")
    draw()
    show()

    wtime2 = clock ( )

    print('')
    print ('    Elapsed wall clock time = %g seconds.'  % ( wtime2 - wtime1 ))
    print ('')
    print ('MD_3D')
    print ('  Normal end of execution.')
    print ('')

    timestamp ( )