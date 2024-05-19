clear
clc
t = 0:0.005:1

R1 = 50E3
R2 = 10E3
R3 = 100E3
C1 = 10E-6
C2 = 50E-9
ei = 1

A = [0 1; -1/(R2*R3*C1*C2) -(1/R1 + 1/R2 + 1/R3)*1/C1]
B = [0; -1/(R1*R2*C1*C2)]
C = [1 0]
D = [0]

ss_system = syslin('c', A, B, C, D)


[yi, xi] = csim('step', t, ss_system)


figure(1)
plot    (t,xi(1, :))
legend  ("$v_o$")

figure(2)
tf_system = ss2tf(ss_system)
evans(tf_system)

figure(3)
bode(ss_system,0.1,100)
