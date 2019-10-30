%% P. Acosta & J. Timana
clc
clear

%% Casos
errorMax=1e-6;
N=200;

%Las estructuras 'Punto1.mat', 'Punto2.mat', 'Punto3.mat' y
%'MarcoTeoric.mat' contiene los datos de el método de N-R como la matriz de
%admitancias 'Ybus' el flag start 'x0', la cantidad de nodos de generación
%'nodosG', la cantidad de nodos de carga 'nodosL', el vector con los
%voltajes de generación 'V0', las cargas reactivas en los nodos de
%generación 'QLoadGen', y el vector solución 'b'.

load('Punto1.mat') %Punto1
[Resul1]=NewtonRaphson(Punto1.Y,Punto1.x0,Punto1.nodosG,Punto1.nodosL,N,errorMax, Punto1.V0, Punto1.QLoadGen,Punto1.b);

load('Punto2.mat') %Punto2
[Resul2]=NewtonRaphson(Punto2.Y,Punto2.x0,Punto2.nodosG,Punto2.nodosL,N,errorMax, Punto2.V0, Punto2.QLoadGen,Punto2.b);

load('Punto3.mat') %Punto3
[Resul3]=NewtonRaphson(Punto3.Y,Punto3.x0,Punto3.nodosG,Punto3.nodosL,N,errorMax, Punto3.V0, Punto3.QLoadGen,Punto3.b);

load('MarcoTeoric.mat') %Marco teórico
[ResulMarco]=NewtonRaphson(MarcoTeoric.Y,MarcoTeoric.x0,MarcoTeoric.nodosG,MarcoTeoric.nodosL,N,errorMax,MarcoTeoric.V0,MarcoTeoric.QLoadGen,MarcoTeoric.b);


%% Flujo de carga

function [s]= NewtonRaphson(Y,x0,nodosG,nodosL,N,errorMax, V0, QLoadGen,b) %Calcula el flujo de carga utilizando el método de Nweton-Raphson

%'s' es una estructura de datos con el tiempo de solución 't', el vector
%solución 'x', las perdidas del sistema 'Perdidas', el vector comparación
%'bt' que es la solución 'b' calculada con 'x0' y el número de iteraciones
%'n'-

%Parametros
%'Y', matriz Ybus; 'x0', vector con solución inicial; 'nodosG', número de
%nodos de generación; 'nodosL', números de nodos de generación; 'N', número
%máximo de iteraciones; 'errorMax', error máximo (i.e criterio de parada),
%'V0' voltajes iniciales de los generadores, 'QLoadGen', es un vector si
%los nodos de generación tienen carga se le resta a Q cálculado para
%mantener las condiciones del flujo de carga.

n=0; %Inicializa el contador en cero
error=1; %Inicializa el error relativo en 1
nodosSlack=1; %Generaliza el nodo 1 como el slack
nodosG=nodosG-1; %NodosG+NodosSlack=Total nodos de generación
nodos=nodosG+nodosL+nodosSlack;%Generaliza los nodos totales
s=struc([]);%Crea la estructura de salida
tic%Inicializa el contador de la solución del algoritmo
while n<N && errorMax<error %Criterios de parada de N-W
    
    %vector de Voltaje de todos los nodos V1 --- VNodos
    V=[V0(1:nodosG+nodosSlack); x0(nodos:end)];
    
    %Thetas vector de desfase respecto al slack de todos los nodos
    %Theta_1------Theta_nodos
    Theta=[0; x0(1:nodos-1)];

    %Vector Potencias de todos los nodos P1---P_nodos
    P=zeros(nodos,1);
    
    %Vector Potencias reactivas de los nodos Q_1-----Q_nodos
    Q=zeros(nodos,1);
    
    %calculo de potencias de todos los nodos
    for i=1:nodos
        for k = 1:nodos
          P(i)= P(i)+abs(V(i))*abs(V(k))*(real(Y(i,k))*cos(Theta(i)-Theta(k)) +imag(Y(i,k))*sin(Theta(i)-Theta(k)));
        end
    end
    
    %calculo de potencias reactancivas de todos los nodos
    for i=1:nodos
        for k=1:nodos 
          Q(i)=Q(i)+ abs(V(i))*abs(V(k))*(real(Y(i,k))*sin(Theta(i)-Theta(k)) - imag(Y(i,k))*cos(Theta(i)-Theta(k)));            
        end 
    end
    
    Q=Q-QLoadGen; %Si los nodos de generación tienen cargas reactivas le resta la carga en el nodo
    %Solución temporal bt=[P_2 ---- P_nodos Q_(NodosG+NodosSlack+1)-----Q_nodos]    
    bt=zeros(nodos,1);
    bt(1:(nodosG+nodosL))=P(2:nodos); %actualiza las potencias del vector de solución tempotal
    bt((nodosG+nodosL+1):nodosG+2*nodosL)=Q((nodosG+2):nodos); %actuazliza las potencias reactivas del vector de solución temporal
    
    %Jacobiano
    J=zeros(nodosG+2*nodosL); %Inicializa el Jacobiano en cero
    MDE=zeros(nodosG+2*nodosL); %Matriz de edición, cada vez que se actuzalice una posición del jacobiano esta matriz indica si se modifico la misma posicion
    
    %Cuando p==q
    %---------------------------------------------------------------------
    %Derivadas parciales de potencias activas {2-nodos} respecto a las
    %derivdas parciales theta {2-nodos}
    
    for i=1:nodosL+nodosG
        J(i,i)=-Q(i+1)- imag(Y(i+1,i+1))*abs(V(i+1))^2;
        MDE(i,i)=MDE(i,i)+1;
    end
    
    %---------------------------------------------------------------------
    %Derivadas parciales de potencias activas {2-nodos} respecto a las
    %derivadas de |V_(nodosGen+nodosSlack+1)| hasta |V_nodos|
    for i=nodosG+1:nodosG+nodosL
        J(i,i+nodosL)=P(i+1)/abs(V(i+1)) + real(Y(i+1,i+1))*abs(V(1+i));
        MDE(i,i+nodosL)=MDE(i,i+nodosL)+1;
    end
    
    %---------------------------------------------------------------------
    %Derivadas parciales de Potencia reactiva {Q_(nodosGen+nodosSlack+1)
    %hasta nodos} respecto a sus derivaadas parciales theta {2-nodos}
    for i = nodosG+nodosL+1:2*nodosL+nodosG
        J(i,i-nodosL)=P(i-nodosL+1)-real(Y(i-nodosL+1,i-nodosL+1))*abs(V(i-nodosL+1))^2;
        MDE(i,i-nodosL)=MDE(i,i-nodosL)+1;
    end
    
    %Derivadas partiales de potencias reactiva Q_(nodosGen+nodosSlack+1)
    %respecto a sus derivadas parciales de |V| desde
    %{(nodosGen+nodosSlack+1) hasta nodos
    
    for i=nodosG+nodosL+1:2*nodosL+nodosG
        J(i,i)=Q(i-nodosL+1)/abs(V(i-nodosL+1)) - imag(Y(i-nodosL+1,i-nodosL+1))*abs(V(i-nodosL+1));
        MDE(i,i)=MDE(i,i)+1;
 
    end
    
    %cuando q!=p
    
    %----------------------------------------------------------------------
    %derivada parcial de potencia respecto a thetha
    for i=1:nodosG+nodosL %recorre 8 ecuaciones de potencia 
        for j=1:nodosG+nodosL %recorre 8 variables Theta. Check
            if i~=j
            J(i,j)=abs(V(i+1)*V(j+1))*(real(Y(i+1,j+1)*sin(Theta(i+1)-Theta(j+1)))-imag(Y(i+1,j+1))*cos(Theta(i+1)-Theta(j+1)));
            MDE(i,j)=MDE(i,j)+1;
            end
        end   
    end
    
    %----------------------------------------------------------------------
    %derivada parcial de potencia reactiva respecto a  Theta
    for i=nodosL+nodosG+1:2*nodosL+nodosG %recorre las 6 ecuaciones de potencia reactiva
        for j=1:nodos-1%recorre las 8 variables Theta
            if i-nodosL~=j
            J(i,j)=-abs(V(i-nodosL+1)*V(j+1))*(real(Y(i-nodosL+1,j+1)*cos(Theta(i-nodosL+1)-Theta(j+1)))-imag(Y(i-nodosL+1,j+1))*sin(Theta(i-nodosL+1)-Theta(j+1)));
            MDE(i,j)=MDE(i,j)+1;
            end
        end
    end
    
    %----------------------------------------------------------------------
    %Derivada parcial de la potencia respecto a los voltajes
    for i=1:nodosL+nodosG %Recorre las 4 ecuacuaciones de potencia
        for j=nodosG+nodosL+1:2*nodosL+nodosG %recorre las 6 variables de voltaje
            if i~=j-nodosL
            J(i,j)=abs(V(i+1))*(real(Y(i+1,j-nodosL+1)*cos(Theta(i+1)-Theta(j-nodosL+1)))+imag(Y(i+1,j-nodosL+1))*sin(Theta(i+1)-Theta(j-nodosL+1)));
            MDE(i,j)=MDE(i,j)+1;  
            end
        end
    end
    
    %---------------------------------------------------------------------
    %Derivada parcial de la potencia reactiva respecto a los voltajes
    for i=nodosG+nodosL+1:2*nodosL+nodosG %Recorre las 6 ecuaciones de potencia reactiva
        for j=nodosG+nodosL+1:2*nodosL+nodosG %Recorre las 6 variables de voltaje
            if j~=i
                J(i,j)=abs(V(i-nodosL+1))*(real(Y(i-nodosL+1,j-nodosL+1))*sin(Theta(i-nodosL+1)-Theta(j-nodosL+1)) -imag(Y(i-nodosL+1,j-nodosL+1))*cos(Theta(i-nodosL+1)-Theta(j-nodosL+1))); 
                MDE(i,j)=MDE(i,j)+1;
            end
        end
        
    end
    
  
    %%%Fin cálculo del Jacobiano%%%
  
    x=x0-J\(bt-b); %Calcula x de la iteración n+1 partiendo de la iteración N.
    error=max(abs(x0-x)); %Calcula el error tomando el valor absoluto infinito de la diferencia de x0 con x
    x0=x; %actualiza el vector 
    n=n+1; %actualiza las iteraciones
end
    t=toc; %tiempo solución método iterativo 
    comp=[bt b]; %Solución: bt calculado respecto a x (solución del método iterativo) y vector solución del problema
    Perdidas=sum(P); %Sumatoria de Pi de {1-nodos}, calcula las perdidas de red.
    A=[(1:nodos)', P, Q, V, 180*Theta/pi];

    %Crea la tabla de resultados
    T = array2table(A,...
    'VariableNames',{'Numero_Nodo','Potencia_Activa_pu','Potencia_Reactiva_pu', 'Voltaje_pu', 'Theta_Voltaje_grados'});
    disp('----------------------------------------------')
    disp('Solución Newton-Raphson')
    disp(T)
    %actualiza la estructura de salida
    s.x=x;
    s.t=t;
    s.Perdidas=Perdidas;
    s.comp=comp;
    s.n=n;
    %%$
    
    % Muestra resultados del método
    disp('El tiempo de ejecución fue de: ')
    disp([num2str(t) ' s'])
    disp('Las perdidas del sistema fueron')
    disp([num2str(Perdidas) ' p.u:'])
    disp(['El sistema se calculó en ' num2str(n) ' iteraciones'])
    disp('----------------------------------------------')
    
end
    
