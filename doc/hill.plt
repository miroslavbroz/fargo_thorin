#!/usr/bin/gnuplot

M_E = 5.97e24
M_S = 1.99e30
au = 1.496e11

a = 1.0*au
a_ = 1.25*au

m = 0.5*M_E

m_ = m
q = m/M_S
q_ = m_/M_S

R_Hill = a*(q/3.)**(1./3.)
R_Hill_ = a_*(q_/3.)**(1./3.)
R_mH = 0.5*(a+a_)*((q+q_)/3.)**(1./3.)

print "a = ", a/au, " au"
print "a_ = ", a_/au, " au"
print "m = ", m/M_E, " M_E"
print "m_ = ", m_/M_E, " M_E"
print "R_H = ", R_Hill/au, " au"
print "R_H' = ", R_Hill_/au, " au"
print "R_mH = ", R_mH/au, " au <-- mutual Hill radius"
print "a'-a = ", (a_-a)/R_mH, " R_mH <-- mutual separation"

multiple = 2.0

a_(a) = a*(1 + multiple*0.5*((q+q_)/3.)**(1./3.))/(1 - multiple*0.5*((q+q_)/3.)**(1./3.))

a1 = 1.0*au
a2 = a_(a1)
a3 = a_(a2)
a4 = a_(a3)
a5 = a_(a4)
a6 = a_(a5)
a7 = a_(a6)
a8 = a_(a7)
a9 = a_(a8)
a10 = a_(a9)

print ""
print "multiple = ", multiple
print "a1 = ", a1/au, " au"
print "a2 = ", a2/au, " au"
#print "a3 = ", a3/au, " au"
#print "a4 = ", a4/au, " au"
#print "a5 = ", a5/au, " au"
#print "a6 = ", a6/au, " au"
#print "a7 = ", a7/au, " au"
#print "a8 = ", a8/au, " au"
#print "a9 = ", a9/au, " au"
#print "a10 = ", a10/au, " au"

#Nrad = 4096
#Nrad = 2048
Nrad = 1024
Nsec = 1536
Rmin = 0.2*au
Rmax = 2.0*au
dR = (Rmax-Rmin)/Nrad
dphi = 2.*pi/Nsec
R = 1.0*au

print ""
print "Nrad = ", Nrad
print "Nsec = ", Nsec
print "Rmin = ", Rmin/au
print "Rmax = ", Rmax/au
print "dR = ", dR/au, " au = ", dR/R_Hill, " R_Hill"
print "R dphi = ", R*dphi/au, " au = ", R*dphi/R_Hill, " R_Hill"

e = 0.001
Q = a*(1+e)
q = a*(1-e)
print "e = ", e
print "Q-q = ", (Q-q)/au, " au = ", (Q-q)/R_Hill, " R_Hill"


