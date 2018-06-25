#!/usr/bin/gnuplot

# units:
# [G] = 1
# [r] = au
# [M] = M_S
# v_kepl = sqrt(GM/r) = 1 @ 1 au
# P = 2pi r/v_kepl = 2pi sqrt(r^3/(GM)) = 2 pi @ 1 au

G = 1
M = 1
r = 1.0
v_kepl = sqrt(G*M/r)
P = 2.*pi*r/v_kepl
dt = P/20
Ntot = 160000
timespan = Ntot*dt

print "code units:"
print "r = ", r
print "P = ", P
print "dt = P/20 = ", dt
print "Ntot = ", Ntot
print "timespan = ", timespan, " = ", timespan/P, " P"

G = 6.67e-11  # SI
M = 1.989e30  # kg
au = 1.496e11  # m
day = 86400.  # s
yr = 365.25*day
kyr = 1.e3*yr
r = r*au
v_kepl = sqrt(G*M/r)
frac = P/dt
P = 2.*pi*r/v_kepl
dt = P/frac
timespan = Ntot*dt

print "\nphysical units:"
print "r = ", r/au, " au"
print "P = ", P/day, " day = ", P/yr, " yr"
print "dt = ", dt/day, " day = ", dt/yr, " yr = ", dt/P, " P"
print "timespan = ", timespan/kyr, " kyr = ", timespan/P, " P"

Ninterm = 500
Nrad = 1024
Nsec = 1536
Nfields = 9
Ngrid = 1

Noutputs = 1.0*Ntot/Ninterm
size1 = 1.0*8*Nrad*Nsec*Nfields
totsize = Ngrid*Noutputs*size1

kB = 1024
MB = 1024*kB
GB = 1024*MB

print "\nfile sizes:"
print "Nrad = ", Nrad
print "Nsec = ", Nsec
print "Ninterm = ", Ninterm
print "Noutputs = ", Noutputs
print "size1 = ", size1, " bytes = ", size1/MB, " MB = ", size1/GB, " GB"
print "totsize = ", totsize, " bytes = ", totsize/MB, " MB = ", totsize/GB, " GB"


