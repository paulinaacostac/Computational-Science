clc
clear 
h=0.000001;%%Separacion de los puntos                                          
t1 = 0:h:30e-3;%%Separacion de los puntos        
%%t1 = 0:h:0.3; 
x = zeros(1,length(t1)); %%Solucion para u_1
y = zeros(1,length(t1)); %%Solucion para u_2

u_s=5*sin(2*pi*50*t1);

y(1) = 0; %% Solucion inicial                                        
x(1) = 0;

R=10000; %Parametros del circuito
C1=0.001;
C2=0.001;
Is=10e-8;
Vt=0.026;


F_xt = @(t,xn,yn)  5*2*pi*50*cos(2*pi*50*t) + (Is/C1)*(exp(-xn/Vt)-1)-(Is/C1)*(exp((xn-yn)/Vt)-1); %Derivadas de funciones (ecuacion diferencial) 

F_yt = @(t,xn,yn)  (exp((xn-yn)/Vt)-1)*(Is/C2)-yn*(1/(R*C2)); 


for i=1:(length(x)-1) %Evaluacion de las funciones con coeficientes de Taylor                              
    k_1 = h*F_xt(t1(i),x(i),y(i));
    l_1 = h*F_yt(t1(i),x(i),y(i));
    
    k_2 = h*F_xt(t1(i)+h/2,x(i)+(k_1)/2,y(i)+(l_1)/2);
    l_2 = h*F_yt(t1(i)+h/2,x(i)+(k_1)/2,y(i)+(l_1)/2);
    
    k_3 = h*F_xt(t1(i)+h/2,x(i)+(k_2)/2,y(i)+(l_2)/2);
    l_3 = h*F_yt(t1(i)+h/2,x(i)+(k_2)/2,y(i)+(l_2)/2);
    
    k_4 = h*F_xt(t1(i)+h,x(i)+k_3,y(i)+l_3);
    l_4 = h*F_yt(t1(i)+h,x(i)+k_3,y(i)+l_3);
    
    x(i+1) = x(i)+(1/6)*(k_1+2*k_2+2*k_3+k_4);
    y(i+1) = y(i)+(1/6)*(l_1+2*l_2+2*l_3+l_4); 
end
figure
plot(t1,x)
grid on;
title('Voltaje en nodo $u_1$','Interpreter','latex');
xlabel('Tiempo (segundos)','Interpreter','latex');
ylabel('$u_1(t)$','Interpreter','latex');



figure
plot(t1,y)
title('y vs T')
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
grid on;
title('Voltaje $u_1$ y $u_2$ en estado estable','Interpreter','latex');
xlabel('Tiempo (segundos)','Interpreter','latex');
ylabel('Voltaje (V)','Interpreter','latex');
legend({'$u_1(t)$','$u_2(t)$'},'Interpreter','latex','location','best');

