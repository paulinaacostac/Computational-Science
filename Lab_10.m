%% J.Timana & P. Acosta 
clc
clear
%% Parámetros
L=150*10^-3; %Inductancia
C=100*10^-6; %Capacitancia
R=12;  %Resistencia
f=50; %frecuencia en Hz

%% Punto 2

syms y(t) %Variable simbolica
Vin=100*sin(2*f*pi*t); %Entrada en terminos simbolicos
ecuacion= y== Vin - (L*C*diff(y,t,2))-(R*C*diff(y,t)); %Ecuacion diferencial circuito RLC
Dy= diff(y,t); %Primera Derivada de y
cond =[y(0)==0, Dy(0)==0]; %Condiciones iniciales
sol(t)= dsolve(ecuacion,cond); %Solución analitica del circuito RLC por solve()


%% Punto 3

yo=[0;0]; %Condiciones iniciales [y(0);Dy(0)]
to=0; %Tiempo inicial 
h=[0.25e-3; 0.5e-3; 1e-3; 2e-3]; %vector h de separación
tmax=[0.4; 0.4; 0.8; 0.2 ]; %vector tiempo


for i=1:length(h)
t=to:h(i):tmax(i); %Tiempo
Yout= sol(t); %Resultado solución analitica
[time, y]=euler2(@A,to,h(i),tmax(i),yo); %Se utiliza para obtener el tiempo y la magnitud
figure
plot(time, y)
xlabel ('Tiempo [s]')
ylabel ('Voltaje [V]')
hold on
plot(t,Yout)
set(gcf,'color','w');
cadena=['Solución con h= ' num2str(h(i)) ' y t_{max}= ' num2str(tmax(i)) 's'];
title(cadena)
legend('Solución método de Euler', 'Solución con función solve()')
end


%%

function [out] = A(t,w)
L=150*10^-3; %Valor de inductancia
C=100*10^-6; %Valor de Capacitancia
R=12;  %Valor de resistencia
out = [w(2); (100*sin(100*pi*t)-R*C*w(2)-w(1))/(L*C)]; %Se utiliza la formula del circuito RLC
end


