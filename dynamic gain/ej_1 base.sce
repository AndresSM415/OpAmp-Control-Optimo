clear;clc;
//Datos de simulación
t0      =   0
tf      =   5
paso    =   0.01
ti      =   (t0:paso:tf)'                   //tiempo para simulación

A       =   [0 1; 2 -1]
B       =   [0  ; 1]
H       =   [0]
Q       =   [1 0; 0 0]
R       =   0.005
r_set   =   [1  ; 0]

k0      =   [0;0;0;0;0]
x0      =   [10 ;-10]


//Ec de Riccati
function dk=ecriccati(ti, k)
    K = [k(1) k(2);k(2) k(3)];
    S = [k(4); k(5)];
    K_p =   K*A  + A'*K + Q - K*B*inv(R)*B'*K
    dk(1:3) = [K_p(1,1) K_p(1,2) K_p(2,2)]
    dk(4:5) =   A'*S - K*B*inv(R)*B'*S - Q*r_set
endfunction
//Simulacion de modelo
k   =   (ode('rk', k0, t0, ti, ecriccati))' //solucion de la ec de estado
t   =   tf - ti + t0                        //tiempo real

k   =   k($:-1:1,:)                         // Reordena ganancias
t   =   t($:-1:1,:)                         // Reordena tiempo
K   =   [k(1,1) k(1,2); k(1,2) k(1,3)]
S   =   [k(1,4)       ; k(1,5)]

//Modelo de estado
function dx =   planta(t, x)                //Ecs. de estado
    x_ss    =   [x(1); x(2)]
    u       =   -inv(R)*B'*K*x - inv(R)*B'*S
    dx(1:2) =   A*x + B*u
endfunction
//Simulacion de modelo
x   =   (ode('rk',x0,t0,t,planta))'         //solucion de la ec de estado
/*
//entrada y funcion de costo
u   =   [] 
dJ  =   []
J   =   zeros()
for i       =   t0+1 : tf/paso+1
    x_opt   =   [x(i,1); x(i,2)]
    xr_opt  =   x_opt - r_set
    u (i)   =   -inv(R)*B'*K*x_opt - inv(R)*B'*S
    dJ(i)   =   (xr_opt'*Q*xr_opt  + R*u(i)^2)/2
    J (i+1) =   J(i) + dJ(i)*paso
end
*/
//Gráficas
figure(1)
plot    (t,k(:,1),      t,k(:,2),   t,k(:,3),   t,k(:,4),   t,k(:,5)    )
legend  ('$K_{11}$',    '$K_{12}$', '$K_{22}$', '$S_{1}$',  '$S_{2}$'   ,3)
xtitle  ('Solución a la Ec. de Riccati',"Tiempo t(s)",'$k(t), s(t)$'    )
xgrid

figure  (2)
N   =   size(t,1)
plot    (t,x(:,1),  t,x(:,2),   t,ones(N,1), '--'   )
legend  ('$x_1$',   '$x_2$',    '$r_1(t)$'          )
xtitle  ('Estados del sistema','Tiempo t (s)','x(t)')
xgrid
/*
figure  (3)
plot    (t,u                                        )
xtitle  ('Entrada del sistema','Tiempo t (s)','u(t)')
xgrid

figure  (4)
plot    (t, J(1:$-1, 1)'                                )
xtitle  ("Funcional de costo", "$t(s)$", "$J(t,x,u)$"   )
xgrid

*/
