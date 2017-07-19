import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import scipy.constants as c

#Question 3.a.i
def energy_n(n):
    m = 9.10938356 * 10**-31
    h = 1.054571800 * 10**-34
    e = 1.6021766208 * 10**-19
    e0 = 8.854187817 * 10**-12
    part1 = m / (2*h*h)
    part2 = (e**2) / (4*math.pi*e0)
    part3 = 1.0 / (n**2.0)
    return -1*part1*part2*part2*part3 / e

#Question 3.a.ii                             
def degToRad(deg):
    return deg*(math.pi/180.0)

def radToDeg(rad):
    return rad*(180.0/math.pi)

#Question 3.a.iii
def sphericalToCartesian(r,theta,phi):
    x = r * math.sin(phi) * math.cos(theta)
    y = r * math.sin(phi) * math.sin(theta)
    z = r * math.cos(phi)
    return x, y, z

def cartesianToSpherical(x, y, z):
    r=np.sqrt(x**2+y**2+z**2)
    rxy=np.sqrt(x**2+y**2)
    if x!=0:
        phi=np.arctan(y/x)
    else:
        phi=math.pi/2
        
    if z!=0 and r!=0:
        theta=np.arccos(z/r)
    else:
        theta=0
    return r, theta, phi

#Question 3.b.i
def fact(n):
    if n == 0:
        return 1.0
    elif (n < 0):
        print ("Input positive/zero number as n")
        return
    else:
        return float(n * fact(n-1))

#Question 3.b.ii
def p00(theta):
    return 1

def p01(theta):
    return np.cos(theta)

def p02(theta):
    return 0.5*(3*np.cos(theta)**2-1)

def p03(theta):
    return 0.5*(5*np.cos(theta)**3-3*np.cos(theta))

def p11(theta):
    return np.sin(theta)

def p12(theta):
    return 3*np.sin(theta)*np.cos(theta)

def p13(theta):
    return 1.5*np.sin(theta)*(5*np.cos(theta)**2-1)

def p22(theta):
    return 3*np.sin(theta)**2

def p23(theta):
    return 15*np.sin(theta)**2*np.cos(theta)

def p33(theta):
    return 15*np.sin(theta)*(1-np.cos(theta)**2)


def assocLegendre(m,l):
    if m==0 and l==0:
        return p00
    elif m==0 and l==2:
        return p02
    elif m==1 and l==1:
        return p11
    elif m==3 and l==3:
        return p33
    elif m==0 and l==1:
        return p01
    elif m==2 and l==3:
        return p23
    elif m==2 and l==2:
        return p22
    elif m==1 and l==3:
        return p13
    elif m==1 and l==2:
        return p12
    elif m==0 and l==3:
        return p03
    else:
        return None


#Question 3.b.iii
def l00(x):
    return 1

def l01(x):
    return -x+1

def l11(x):
    return -2*x+4

def l10(x):
    return 1

def l02(x):
    return x*x-4*x+2

def l12(x):
    return 3*x*x-18*x+18

def l22(x):
    return 12*x*x-96*x+144

def l21(x):
    return -6*x+18

def l02(x):
    return 2

def l03(x):
    return -x*x*x+9*x*x-18*x+6

def l13(x):
    return -4*x*x*x+48*x*x-144*x+96

def l32(x):
    return 60*x*x-600*x+1200

def l33(x):
    return -120*x*x*x+2160*x*x-10800*x+14400

def l23(x):
    return -20*x*x*x+300*x*x-1200*x+1200

def l31(x):
    return -24*x+96

def l30(x):
    return 6

def assocLaguerre(p,qmp):
    if (p == 0):
        if (qmp == 0):
            return l00
        elif (qmp == 1):
            return l01
        elif (qmp == 2):
            return l02
        elif (qmp == 3):
            return l03

    elif (p == 1):
        if (qmp == 0):
            return l10
        elif (qmp == 1):
            return l11
        elif (qmp == 2):
            return l12
        elif (qmp == 3):
            return l13

    elif (p == 2):
        if (qmp == 0):
            return l20       
        elif (qmp == 1):
            return l21
        elif (qmp == 2):
            return l22
        elif (qmp == 3):
            return l23

    elif (p == 3):
        if (qmp == 0):
            return l30       
        elif (qmp == 1):
            return l31
        elif (qmp == 2):
            return l32
        elif (qmp == 3):
            return l33

    else:
        return ("Non-valid input")

#Question 3.c.i    
def azimuth(m,phi):
	return np.exp(1j* m * phi)

def angular_wave_func(m,l,theta,phi):
	if m>=0:
		eps=(-1)**m
	else:
		eps=1
	y1=(2*l+1)/(4.0*c.pi)
	y2=float(fact(l-abs(m)))/fact(l+abs(m))
	ymid=np.sqrt(y1*y2)
	pfunc=assocLegendre(m,l)
	x=azimuth(m,phi)
	#print ymid,x, pfunc(theta)
	return eps*ymid*x*pfunc(theta)

#Question 3.c.ii
def radial_wave_func(n,l,r):
    a = 5.29177*(10**(-11))
    part1a = (2.0/(n*a))**3.0
    part1b = (fact(n-l-1)) / (2*n*((fact(n+l))**3.0))
    part1 = (part1a*part1b)**0.5
    part2 = np.exp(-r/(n*a))
    part3 = ((2*r)/(n*a))**l
    part4 = assocLaguerre((2*l+1),(n-l-1))
    partfinal = part1*part2*part3*part4((2*r)/(n*a))
    return partfinal

#Question 3.c.iii
#For the density plots we use (r^2)(Rnl^2)
a = 5.29177*(10**(-11))
x = np.linspace(1*10**-11,1.2*10**-9,200)
y=[]

for b in x:
    y.append(radial_wave_func(3,1,b) / a**-1.5 )

for b in range(0,len(x)-1):
    x[b] = x[b] / a

plt.figure(1)
plt.title("Radial Wave Function for n = 3 and l = 1")
plt.xlabel("Distance From Radius")
plt.ylabel("Probabilty")
plt.axis('auto')
plt.plot(x,y,'ro')
plt.show()

###Question 3.d.i
def hydrogen_wave_func(n, m, l, roa, Nx, Ny, Nz):
	x=np.linspace(-roa,roa,Nx)
	y=np.linspace(-roa,roa,Ny)
	z=np.linspace(-roa,roa,Nz)
	xx,yy,zz=np.meshgrid(x,y,z)
	#print len(xx), len(yy),len(zz)
	vfunc=np.vectorize(cartesianToSpherical)
	roaa, theta, phi=vfunc(xx,yy,zz)
	a=c.physical_constants['Bohr radius'][0]
	r=roaa*a
	rmag=radial_wave_func(n,l,r)
	angumag=angular_wave_func(m,l,theta,phi)
	final=np.absolute(rmag*angumag)**2
	return xx,yy,zz,final
    
##x = np.linspace(-1.2*10**-9,1.2*10**-9,10)
##y = np.linspace(-1.2*10**-9,1.2*10**-9,10)
##z = np.linspace(-1.2*10**-9,1.2*10**-9,10)
##
##xnew = []
##ynew = []
##znew = []
##
##result = []
##
##for a in x:
##    for b in y:
##        for c in z:
##            r, theta, phi = cartesianToSpherical(a,b,c)
##            xnew.append(a)
##            ynew.append(b)
##            znew.append(c)
##            part1 = cmath.polar(((angular_wave_func(0,0,theta,phi))*(radial_wave_func(2,0,r)))**2)
##            result.append(part1[0])
##
##
##for a in range(0,len(result)):
##    result[a] = result[a] / max(result)
##
    
x,y,z,mag=hydrogen_wave_func(2,0,1,10,20,20,20)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for a in range(0,len(mag)):
    for b in range(0,len(mag)):
        for c in range(0,len(mag)):
            ax.scatter(x[a][b][c],y[a][b][c],z[a][b][c], marker='o',alpha=(mag[a][b][c]/np.amax(mag)))

plt.show()
