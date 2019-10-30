clc
clear 
h=0.000001;  %%Separacion de los puntos                                           
t1 = 0:h:0.025; %%Vector de tiempo que determina el numero de puntos
%%t1 = 0.2:h:1;  
x = zeros(1,length(t1)); %%Solucion para u_1
y = zeros(1,length(t1)); %%Solucion para i_L
z = zeros(1,length(t1)); %%Solucion para u_2

x(1)=0;
y(1)=0;
z(1)=0;
u_s=10*sin(2*pi*50*t1); %%Funcion de voltaje de entrada

R=50; %%parametros del circuito
C1=1000e-6;
C2=1000e-6;
L=0.1;
Is=10e-8;
Vt=0.026;

Fx = @(t,xn,yn,zn) (1/C1)*yn - (1/(C1*R))*xn; %Derivadas de las funciones (ecuacion diferencial)
Fy = @(t,xn,yn,zn) (1/L)*zn - (1/L)*xn;
Fz = @(t,xn,yn,zn) (1/C2)*Is*(exp((abs(10*sin(2*pi*50*t))-zn)/(2*Vt))-1)-(1/C2)*yn;

for i=1:length(t1)-1
    %Evaluacion de las funciones con coeficientes de Taylor   
    k1=h*Fx(t1(i),x(i),y(i),z(i));
    l1=h*Fy(t1(i),x(i),y(i),z(i));
    m1=h*Fz(t1(i),x(i),y(i),z(i));
    
    k2=h*Fx(t1(i)+h/2,x(i)+k1/2,y(i)+l1/2,z(i)+m1/2);
    l2=h*Fy(t1(i)+h/2,x(i)+k1/2,y(i)+l1/2,z(i)+m1/2);
    m2=h*Fz(t1(i)+h/2,x(i)+k1/2,y(i)+l1/2,z(i)+m1/2);
    
    k3=h*Fx(t1(i)+h/2,x(i)+k2/2,y(i)+l2/2,z(i)+m2/2);
    l3=h*Fy(t1(i)+h/2,x(i)+k2/2,y(i)+l2/2,z(i)+m2/2);
    m3=h*Fz(t1(i)+h/2,x(i)+k2/2,y(i)+l2/2,z(i)+m2/2);
    
    k4=h*Fx(t1(i)+h,x(i)+k3,y(i)+l3,z(i)+m3);
    l4=h*Fy(t1(i)+h,x(i)+k3,y(i)+l3,z(i)+m3);
    m4=h*Fz(t1(i)+h,x(i)+k3,y(i)+l3,z(i)+m3);
    
    x(i+1)= x(i)+ (1/6)*(k1+2*k2+2*k3+k4);
    y(i+1)= y(i)+ (1/6)*(l1+2*l2+2*l3+l4);
    z(i+1)= z(i)+ (1/6)*(m1+2*m2+2*m3+m4);
    

    
    
end

figure
plot(t1,x)
title('x vs T')
grid on;
title('Voltaje en nodo $u_1$','Interpreter','latex');
xlabel('Tiempo (segundos)','Interpreter','latex');
ylabel('$u_1(t)$','Interpreter','latex');


figure
plot(t1,y)
title('y vs T')
grid on;
title('Corriente en inductancia L $i_L$','Interpreter','latex');
xlabel('Tiempo (segundos)','Interpreter','latex');
ylabel('$i_(t)$','Interpreter','latex');


figure
plot(t1,z)
title('z vs T')
grid on;
title('Voltaje en nodo $u_2$','Interpreter','latex');
xlabel('Tiempo (segundos)','Interpreter','latex');
ylabel('$u_2(t)$','Interpreter','latex');

figure
plot(t1,u_s)
title('U_s')
grid on;
title('Voltaje de entrada','Interpreter','latex');
xlabel('Tiempo (segundos)','Interpreter','latex');
ylabel('$u_s(t)$','Interpreter','latex');

figure
plot(t1,x)
hold on
plot(t1,y)
hold on
plot(t1,z)
grid on;
title('$u_1(t)$, $u_2(t)$, $i_L(t)$ en estado estable','Interpreter','latex');
xlabel('Tiempo (segundos)','Interpreter','latex');
ylabel('Voltaje (V)','Interpreter','latex');
legend({'$u_1(t)$','$i_L(t)$','$u_2(t)$'},'Interpreter','latex','location','best');

