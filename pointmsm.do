mata

function tabc(Y1A1L1, Y1A0L1, Y1A1L0, Y1A0L0, Y0A1L1, Y0A0L1, Y0A1L0, Y0A0L0, w, s, n){

dAtAL1 = (Y1A1L1, Y0A1L1, Y1A0L1, Y0A0L1)'
dataL0 = (Y1A1L0, Y0A1L0, Y1A0L0, Y0A0L0)'
N = (dAtAL1\ dataL0)
L = J(4,1,1)\ J(4,1,0)
A = J(2,1,1)\ J(2,1,0)\J(2,1,1)\ J(2,1,0)
Y = (1,0,1,0,1,0,1,0)'

prA1L1 = (Y1A1L1 + Y0A1L1)/(Y1A1L1 + Y0A1L1 + Y1A0L1 + Y0A0L1)
prA0L1 = (Y1A0L1 + Y0A0L1)/(Y1A1L1 + Y0A1L1 + Y1A0L1 + Y0A0L1)
prA1L0 = (Y1A1L0 + Y0A1L0)/(Y1A1L0 + Y0A1L0 + Y1A0L0 + Y0A0L0)
prA0L0 = (Y1A0L0 + Y0A0L0)/(Y1A1L0 + Y0A1L0 + Y1A0L0 + Y0A0L0)

sumall = (Y1A1L1+Y1A0L1+Y1A1L0+Y1A0L0+Y0A1L1+Y0A0L1+Y0A1L0+Y0A0L0)
prA = (Y1A1L1 + Y0A1L1 + Y1A0L0 + Y0A0L0)/ sumall

prA1L1a = prA
prA0L1a = prA
prA1L0a = prA
prA0L0a = prA

pr = (prA1L1,prA1L1, prA0L1,prA0L1, prA1L0,prA1L0, prA0L0,prA0L0)'
pr2 = (prA1L1a,prA1L1a, prA0L1a,prA0L1a, prA1L0a,prA1L0a, prA0L0a,prA0L0a)'



if (w == 0){
w = J(8,1,1)
N = N
} 
else 
{
if (n == 0 & s == 0){
w = J(8,1,1):/pr
N = N:*w
}

if (n == 1 & s == 0){
w0 = J(8,1,1):/pr
w = (w0 :* sumall) :/ J(8,1,sum(w0 :* N))
N = N:*w
}

if (n == 0 & s == 1){
w = pr2 :/ pr
N = N:*w
}

if (n == 1 & s == 1){
sw = pr2 :/ pr
w = (sw :* sumall) :/ J(8,1,sum(sw :* N))
N = N:*w
}
}

data = (L,A,Y,N, w)
return(data)
}


function causal(Y1A1L1, Y1A0L1, Y1A1L0, Y1A0L0, Y0A1L1, Y0A0L1, Y0A1L0, Y0A0L0, w, s, n){

dAtAL1 = (Y1A1L1, Y0A1L1, Y1A0L1, Y0A0L1)'
dataL0 = (Y1A1L0, Y0A1L0, Y1A0L0, Y0A0L0)'
N = (dAtAL1\ dataL0)
L = J(4,1,1)\ J(4,1,0)
A = J(2,1,1)\ J(2,1,0)\J(2,1,1)\ J(2,1,0)
Y = (1,0,1,0,1,0,1,0)'

prA1L1 = (Y1A1L1 + Y0A1L1)/(Y1A1L1 + Y0A1L1 + Y1A0L1 + Y0A0L1)
prA0L1 = (Y1A0L1 + Y0A0L1)/(Y1A1L1 + Y0A1L1 + Y1A0L1 + Y0A0L1)
prA1L0 = (Y1A1L0 + Y0A1L0)/(Y1A1L0 + Y0A1L0 + Y1A0L0 + Y0A0L0)
prA0L0 = (Y1A0L0 + Y0A0L0)/(Y1A1L0 + Y0A1L0 + Y1A0L0 + Y0A0L0)

sumall = (Y1A1L1+Y1A0L1+Y1A1L0+Y1A0L0+Y0A1L1+Y0A0L1+Y0A1L0+Y0A0L0)
prA = (Y1A1L1 + Y0A1L1 + Y1A0L0 + Y0A0L0)/ sumall

prA1L1a = prA
prA0L1a = prA
prA1L0a = prA
prA0L0a = prA

pr = (prA1L1,prA1L1, prA0L1,prA0L1, prA1L0,prA1L0, prA0L0,prA0L0)'
pr2 = (prA1L1a,prA1L1a, prA0L1a,prA0L1a, prA1L0a,prA1L0a, prA0L0a,prA0L0a)'



if (w == 0){
w = J(8,1,1)
N = N
} 
else 
{
if (n == 0 & s == 0){
w = J(8,1,1):/pr
N = N:*w
}

if (n == 1 & s == 0){
w0 = J(8,1,1):/pr
w = (w0 :* sumall) :/ J(8,1,sum(w0 :* N))
N = N:*w
}

if (n == 0 & s == 1){
w = pr2 :/ pr
N = N:*w
}

if (n == 1 & s == 1){
sw = pr2 :/ pr
w = (sw :* sumall) :/ J(8,1,sum(sw :* N))
N = N:*w
}
}

data = (L,A,Y,N, w)

mL1 = (N[1], N[3])\(N[2], N[4])
mL0 = (N[5], N[7])\(N[6], N[8])
m = (mL1 + mL0)
m = (mL1 + mL0)\colsum(m)

r1 = m[1,1]/m[3,1]
r0 = m[1,2]/m[3,2]
rd = r1 - r0
rr = r1 / r0
or = (r1/(1-r1)) / (r0/(1-r0))
est = (rd,rr,or)

return(est)
}

end