clear;clc;
//Datos de simulación
t0          =   0
tf          =   1
paso        =   0.1
ti          =   (t0:paso:tf)'                   //tiempo para simulación

R1 = 50E3
R2 = 10E3
R3 = 100E3
C1 = 10E-6
C2 = 50E-9

A = [0 1; -1/(R2*R3*C1*C2) -(1/R1 + 1/R2 + 1/R3)*1/C1]
B = [0; -1/(R1*R2*C1*C2)]

H       =   [0.01 0; 0 0]
Q       =   [0.01 0; 0 0.01]
R       =   0.01
r_set   =   [10; 0]

k0      =   [H(1,1); 0; H(2,2); -H*r_set]
x0      =   [0; 0]

//Ec de Riccati
function dk=ecriccati(ti, k)
    //r_set   =   [0.2*(tf-ti+t0); 0]
    K       =   [k(1) k(2); k(2) k(3)];
    S       =   [k(4); k(5)];
    K_p     =   K*A  + A'*K + Q - K*B*inv(R)*B'*K
    dk(1:3) =   [K_p(1,1) K_p(1,2) K_p(2,2)]
    dk(4:5) =   A'*S - K*B*inv(R)*B'*S - Q*r_set
endfunction
//Simulacion de modelo
tic();
k   =   (ode('rk', k0, t0, ti, ecriccati))' //solucion de la ec de estado
disp(toc()); disp(1);

t   =   tf - ti + t0                        //tiempo real
k   =   k($:-1:1,:)                         // Reordena ganancias
t   =   t($:-1:1,:)                         // Reordena tiempo

//Modelo de estado
function dx =   planta(t, x)                //Ecs. de estado
    /*k12 =   interp1(ti,k(:,2),t)
    k22 =   interp1(ti,k(:,3),t)
    s2  =   interp1(ti,k(:,5),t)
    //u       =   -inv(R)*B'*K*x - inv(R)*B'*S
    u   =   -inv(R)*[k12*x(1) + k22*x(2) - s2]
    dx(1:2) =   A*x + B*u*/
    K       =   [0 0; interp1(ti,k(:,2),t) interp1(ti,k(:,3),t)]
    S       =   [0; interp1(ti,k(:,5),t)]
    u       =   -inv(R)*B'*K*x - inv(R)*B'*S
    dx(1:2) =   A*x + B*u
endfunction
tic();
x   =   (ode('rk',x0,t0,t,planta))'         //solucion de la ec de estado
disp(toc()); disp(2);

tic();
//entrada y funcional de costo
u_opt   =   []
J_opt   =   zeros()
for i           =   t0+1 : tf/paso+1
    x_opt       =   [x(i,1); x(i,2)]
    //r_set       =   [0.2*i; 0]
    xr_opt      =   x_opt - r_set
    K           =   [k(i,1) k(i,2); k(i,2) k(i,3)]
    S           =   [k(i,4); k(i,5)]
    u_opt (i)   =   -inv(R)*B'*K*x_opt - inv(R)*B'*S
    dJ_opt      =   (xr_opt'*Q*xr_opt  + R*u_opt(i)^2)/2
    J_opt (i+1) =   J_opt(i) + dJ_opt*paso
end
disp(toc()); disp(3);
//Gráficas
figure(1)
plot    (t,[k(:,1) k(:,2) k(:,3) k(:,4) k(:,5)]    )
legend  ('$K_{11}$', '$K_{12}$', '$K_{22}$', '$S_{1}$', '$S_{2}$' ,3)
xtitle  ('Solución a la Ec. de Riccati',"Tiempo t(s)",'$k(t), s(t)$')
xgrid

figure  (2)
N   =   size(t,1)
plot    (t,x(:,1),  t,x(:,2),  t,r_set(1)*ones(t), '--'   )
legend  ('$x_1$',   '$x_2$',   '$r_1(t)$'          ,2)
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


  /*k11 = k(1); k12 = k(2); k22 = k(3)
    s1 = k(4); s2 = k(5)
    dk(1)   =   -4E3*k12 + Q(1,1) + 16E6*inv(R)*k12^2
    dk(2)   =   -2E3*k22 - 13*k12 + k11 + 16E6*inv(R)*k12*k22
    dk(3)   =   -26*k22 + 2*k12 + Q(2,2) + 16E6*inv(R)*k22^2
    dk(4)   =   -[2000 + 16E6*inv(R)*k12]*s1 - Q(1,1)*r_set(1)
    dk(5)   =   -[13 + 16E6*inv(R)*k22]*s2 + s1 - Q(2,2)*r_set(2)*/
