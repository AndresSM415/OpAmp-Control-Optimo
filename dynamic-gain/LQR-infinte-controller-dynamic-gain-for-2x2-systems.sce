clear;clc;
//Datos de simulación
t0          =   0
tf          =   0.08
paso        =   0.001
ti          =   (t0:paso:tf)'                   //tiempo para simulación

R1 = 50E3
R2 = 10E3
R3 = 100E3
C1 = 10E-6
C2 = 50E-9

A = [0 1; -1/(R2*R3*C1*C2) -(1/R1 + 1/R2 + 1/R3)*1/C1]
B = [0; -1/(R1*R2*C1*C2)]
C = [1 0]
D = [0]
H           =   [0]
Q           =   [1.9 0; 0 0]
R           =   0.1
r_setpoint  =   [10  ; 0]

k0          =   [0; 0; 0; 0; 0]
x0          =   [0 ;0]

//Ec de Riccati
function dk =   ecriccati(ti, k)
    K       =   [k(1) k(2); k(2) k(3)];
    S       =   [k(4); k(5)];
    K_der   =   K*A + A'*K + Q - K*B*inv(R)*B'*K
    dk(1:3) =   [K_der(1,1) K_der(1,2) K_der(2,2)]
    dk(4:5) =   A'*S - K*B*inv(R)*B'*S - Q*r_setpoint
endfunction
k   =   (ode('rk', k0, t0, ti, ecriccati))' //solucion de la ec de estado
k   =   k($:-1:1,:)                         // Reordena ganancias
K   =   [k(1,1) k(1,2); k(1,2) k(1,3)]
S   =   [k(1,4); k(1,5)]
t   =   tf - ti + t0                        //tiempo real
t   =   t($:-1:1,:)                         // Reordena tiempo

//Modelo de estado
function dx =   planta(t, x)                //Ecs. de estado
    u       =   -inv(R)*B'*K*x - inv(R)*B'*S
    dx(1:2) =   A*x + B*u
endfunction
x   =   (ode('rk', x0, t0, t, planta))'     //solucion de la ec de estado

//entrada y funcional de costo
u_opt   =   []
J_opt   =   zeros()
for i           =   t0+1 : tf/paso+1
    x_opt       =   [x(i,1); x(i,2)]
    xr_opt      =   x_opt - r_setpoint
    u_opt (i)   =   -inv(R)*B'*K*x_opt - inv(R)*B'*S
    dJ_opt      =   (xr_opt'*Q*xr_opt  + R*u_opt(i)^2)/2
    J_opt (i+1) =   J_opt(i) + dJ_opt*paso
end

//Gráficas
figure(1)
plot    (t,k(:,1),      t,k(:,2),   t,k(:,3),   t,k(:,4),   t,k(:,5)    )
legend  ('$K_{11}$',    '$K_{12}$', '$K_{22}$', '$S_{1}$',  '$S_{2}$', 3)
xtitle  ('Solución a la Ec. de Riccati',"Tiempo t(s)",'$k(t), s(t)$'    )
xgrid

figure  (2)
plot    (t,x(:,1),  t,x(:,2),   t,r_setpoint(1), '--'   )
legend  ('$x_1$',   '$x_2$',    '$r_1(t)$'          , 4)
xtitle  ('Estados del sistema','Tiempo t (s)','x(t)')
xgrid

figure  (3)
plot    (t,u_opt                                    )
xtitle  ('Entrada del sistema','Tiempo t (s)','u(t)')
xgrid

figure  (4)
plot    (t, J_opt(1:$-1)                                )
xtitle  ("Funcional de costo", "$t(s)$", "$J(t,x,u)$"   )
xgrid
