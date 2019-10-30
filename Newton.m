clc
clear
%%

N=1000;
n=0;
errorMax=1e-6;
error=1;

%   1    2    3   4   5   
Y=[ 
   2-1i*20,   -1+10*1i,     0,         -1+1i*10,   0;
   -1+1i*10,  3-1i*30,      -1+1i*10,   -1+1i*10,   0;
   0,         -1+1i*10,     2-1i*20,   0,          -1+1i*10;
   -1+1i*10,  -1+1i*10,     0,         3-1i*30,    -1+1i*10;
   0,         0,            -1+1i*10,  -1+1i*10,   2-1i*20
   ];
%x0=[Thetag2 Thetag3 Theta4 Theta5 V4 V5]' 
% x0=[0; 0; 0; 0; 1; 1]; %Condiciones iniciales guia
x0=[-5*pi/180; -10*pi/180; -10*pi/180; -15*pi/180; 1; 1]; %Solucion esperada

%Incognitas P1, Q1, Qg2, Thetag2, Qg3, Thetag3, V4, V5,  Theta4, Theta5
%x=[Thetag2 Thetag3 Theta4 Theta5 V4 V5]

%b vector solución b=[P2 P3 P4 P5 Q4 Q5]'
b=[0.8830;(0.2076-0.2); 1.7137; 1.7355; -(0.5983-1); -(0.5496-0.8)];


tic
while n<N && errorMax<error
    
    %vector norma voltajes V1 V2 V3 V4 V5
    V=[1; 1; 1; x0(5:6)];

    %Thetas
    Theta=[0; x0(1:4)];

    %Solución temporal bt=[P2 P3 P4 P5 Q4 Q5]
    bt=zeros(6,1);
    
    
    
    %calculo de potencias
    for i=1:4
        for k = 1:5
            bt(i)=bt(i)+abs(V(i+1))*abs(V(k))*(real(Y(i+1,k))*cos(Theta(i+1)-Theta(k)) +  imag(Y(i+1,k))*sin(Theta(i+1)-Theta(k)));
        end
    end
    
    %calculo de reactancias
    for i=5:6
        for k=1:5   
            bt(i)=bt(i)+abs(V(i-1))*abs(V(k))*(real(Y(i-1,k))*sin(Theta(i-1)-Theta(k)) - imag(Y(i-1,k))*cos(Theta(i-1)-Theta(k)));
        end 
    end
    %Jacobiano
    J=zeros(6);
    MDE=zeros(6);
    %Cuando p==q
    %---------------------------------------------------------------------
    %Derivadas parciales de potencias i en {2-5}  a sus derivadas
    %parciales Theta {2-5} respectivamente. CHeck
    for i=1:4
        J(i,i)= -bt(i)*tan(Theta(i+1))-imag(Y(i+1,i+1))*abs(V(i+1))^2;
        MDE(i,i)=MDE(i,i)+1;
        if i==2
            J(i,i)=J(i,i);
        end
    end
    %---------------------------------------------------------------------
    %Derivadas parciales de potencias i en {4-5} a sus derrivadas parciales
    %Voltaje 4-5 respectivamente. Chech
    for i=3:4
        J(i,i+2)=bt(i)/abs(V(i+1)) + real(Y(i+1,i+1))*V(i+1);
        MDE(i,i+2)=MDE(i,i+2)+1;
    end
    %---------------------------------------------------------------------
    %Derivadas parciales de Potencia reactiva i en {4-5} a sus derivadas
    %parciales Theta {4-5} respectivamente. Chech
    for i = 5:6
        J(i,i-2)=bt(i-2)-real(Y(i-1,i-1))*abs(V(i-1))^2;
        MDE(i,i-2)=MDE(i,i-2)+1;
    end
    %----------------------------------------------------------------------
    %Derivadas partiales de potencias reactivas i en {4-5} a sus dereivadas
    %parciales voltaje {4-5}. Chech
    
    for i=5:6
        J(i,i)=bt(i-2)*tan(Theta(i-1))/abs(V(i-1)) -imag(Y(i-1,i-1))*abs(V(i-1));
        MDE(i,i)=MDE(i,i)+1;
    end
    
    %cuando q!=p
    %----------------------------------------------------------------------
    %derivada parcial de potencia respecto a thetha
    for i=1:4 %recorre 4 ecuaciones de potencia 
        for j=1:4 %recorre 4 variables Theta. Check
            if i~=j
            J(i,j)=abs(V(i+1)*V(j+1))*(real(Y(i+1,j+1)*sin(Theta(i+1)-Theta(j+1)))-imag(Y(i+1,j+1))*cos(Theta(i+1)-Theta(j+1)));
            MDE(i,j)=MDE(i,j)+1;
            end
        end   
    end
    
    %----------------------------------------------------------------------
    %derivada parcial de potencia reactiva respecto a  Theta
    for i=5:6 %recorre las 2 ecuaciones de potencia reactiva
        for j=1:4%recorre las 4 variables Theta
            if i-2~=j
            J(i,j)=-abs(V(i-1)*V(j+1))*(real(Y(i-1,j+1)*cos(Theta(i-1)-Theta(j+1)))-imag(Y(i-1,j+1))*sin(Theta(i-1)-Theta(j+1)));
            MDE(i,j)=MDE(i,j)+1;
            end
        end
    end
    
    %----------------------------------------------------------------------
    %Derivada parcial de la potencia respecto a los voltajes
    for i=1:4 %Recorre las 4 ecuacuaciones de potencia
        for j=5:6 %recorre las 2 variables de voltaje
            if i~=j-2
            J(i,j)=abs(V(i+1))*(real(Y(i+1,j-1)*cos(Theta(i+1)-Theta(j-1)))+imag(Y(i+1,j-1))*sin(Theta(i+1)-Theta(j-1)));
            MDE(i,j)=MDE(i,j)+1;  
            end
        end
    end
    
    %---------------------------------------------------------------------
    %Derivada parcial de la potencia reactiva respecto a los voltajes
    for i=5:6
        for j=5:6
            if j~=i
                J(i,j)=abs(V(j-1))*(real(Y(i-2,i-2))*sin(Theta(i-1)-Theta(j-1)) -imag(Y(i-2,j-2))*cos(Theta(i-1)-Theta(j-1))); 
                MDE(i,j)=MDE(i,j)+1;
            end
        end
        
    end
    
    
    
    x=x0-J\(bt-b);
    error=max(abs(x0-x));
    x0=x;
    n=n+1;
       
end
   tt=toc; 
    
    
    th=[0;x(1:4)];
    v=[1; 1; 1;x(5:6)];
    z=zeros(6,1);
    
    for i=1:4
        for k = 1:5
            z(i)=z(i)+abs(v(i+1))*abs(v(k))*(real(Y(i+1,k))*cos(th(i+1)-th(k)) +  imag(Y(i+1,k))*sin(th(i+1)-th(k)));
        end
    end
    for i=5:6
        for k=1:5   
            z(i)=z(i)+abs(v(i-1))*abs(v(k))*(real(Y(i-1,k))*sin(th(i-1)-th(k)) - imag(Y(i-1,k))*cos(th(i-1)-th(k)));
        end 
    end

    comp=[b z];

xdegree=(180/pi).*Theta
