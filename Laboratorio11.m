%% J. Timana & P. Acosta
clc
clear
R=[]; %Matriz donde se almacena el calculo de Rt con los métodos trapezoidal, Simpson y Gauss-Legendre. Se inicializa como vacia.
n=100; %Grado de aproximación
h=pi/n; %Separación de pasos
u= 0:h:3; %Vector donde se evalua la formula
%%  Trapecio
R=[R;60*Trapecio(u,h,n)]; %Calcula y anade a la matriz el resultado de la integral por el método trapezoidal

%% Simpson
R=[R;60*Simpson(h,n,u)]; %Calcula y anade a la matriz el resultado de la integral por el método de Simpson


%% Gauss Legendre
R=[R;GaussLe(0.00001,pi,n,u)]; %Calcula y anade a la matriz el resultado de la integral por el método de Gauss-Legendre

%% Graficas

%Todas las graficas de los distintos métodos
figure
plot(u,R(1,:),u,R(2,:),u,R(3,:))
set(gcf,'color','w');
xlabel('[u]' )
ylabel('R_r Resistencia de radiación [\Omega]')
legend('Método de trapecio', 'Formula de Simpson', 'Gauss-Legendre')
title('Comparación de los tres métodos')

%Trapecio
figure
plot(u,R(1,:))
set(gcf,'color','w');
xlabel('[u]' )
ylabel('R_r Resistencia de radiación [\Omega]')
title('Método de trapecio')

%Simpson
figure
plot(u,R(2,:))
set(gcf,'color','w');
xlabel('[u]' )
ylabel('R_r Resistencia de radiación [\Omega]')
title('Formula de Simpson')

%Gauss-Legendre
figure
plot(u,R(3,:))
set(gcf,'color','w');
xlabel('[u]' )
ylabel('R_r Resistencia de radiación [\Omega]')
title('Gauss-Legendre')


%% Formula de Simpson
function [ R ] = Simpson(h,n,u)
a=0 ; %Inicializa el límite inferior en 0
sum= 0; %Sumatoria que almacena el valor parcial de la inegral
LimInfe=0.00001; %Manejo de errores

%Calcula la sumatoria
for i=1:(n/2)-1
    sum=sum+Fun(u, a+i*h); %a+i*h se utiliza para calcular el respectivo xj
end
sum2=0;
%Calcula la segunda sumatoria
for i=1:(n/2)
    sum2=Fun(u,a+((2*i)-1)*h)+sum2; 
end

R= (h/(3))*(Fun(u,LimInfe)+2*sum+4*sum2 + Fun(u,pi)); %Resultado de la integral
end


%% Argumento integral
function [Rt] = Fun(u,theta) 
Rt= ((cos(pi*u*cos(theta))-cos(pi*u)).^2)./sin(theta); %Calcula el argumento de la integral
end

%% Método del trapecio
function [R] = Trapecio(u,h,n)
a=0; 
sum=0;
for i=1:n-1
   sum=sum+Fun(u,a+i*h); 
R= (Fun(u,0.00001)+2*sum+ Fun(u,pi))*(h/(2)); %Realiza dentro del orden dado la formula del trapecio
end
end
%%

function [R] = GaussLe(a,b,n,u)
R=0; %Se inicializa el resultado de la integral en 0
syms x %Variable simbólica
k=(legendreP(n,x));%Polinomio de legrenge de grado n
xi=double(solve(k==0)).'; %Solucionar para las raices del polinomio
c=zeros(1,n); % Inicializa el vector C en 0
%Buscar cada Ci
for i=1:n
c(i)=(2*(1-(xi(i).^2))/((n^2)*(legendreP(n-1,xi(i)).^2)));
end
%Calcula la integral de grado n
for i=1:n
   R=R+c(i)*Fun(u,((b-a)*xi(i)+a+b)/2)*(b-a)*0.5*60;
end
end