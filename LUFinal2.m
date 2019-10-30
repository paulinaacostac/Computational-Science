%% P. Acosta & J. Timana: Lab4-Computaci�n
clc
clear
%% Segundo Punto
%MatrizA
A=[1 1 0 3;2 1 -1 1;3 -1 -1 2;-1 2 3 -1];
F=[2.121 -3.460 0 5.217; 0 5.193 -2.197 4.206; 5.132 1.414 3.141 0;-3.111 -1.732 2.718 5.212]; %Matriz A
b=[1; 2; 3]; %Vector columna con resultados aleatorios
[A1, u1]=MatrizLU(A); % Ejecuta el algoritmo de matriz LU. 'A1' matriz LU, 'u1' matriz de pivotes parciales calculados sobre 'A'
[L,U,P]=factorizacionLU(A1,u1);% Ejectuta el algoritmo para obtener las matrices L,U,P. 'L' matriz diagonal inferior, 'U' matriz diagonal superior, 'P' matriz de permutaci�n.
y=L\(P*b); %Se soluciona el sistema Ly=Pb. 'L'\('P'*'b')=inv(L)*P*b. 'y' es el vector soluci�n.
x=U\y;%Se soluciona el sistema Ux=b. 'U'\'y'=inv(U)*y. 'x' es el vector soluci�n.
xt=A\b; %Soluci�n t�orica Ax=y. 'A'\'b'=inv(A)*b. 'xt' es la soluci�n te�rica.
comp1=[x xt]; %Se unen en una matriz los vectores columnas con la soluc�n mediante el algoritmo y la te�rica respectivamente
error1=abs(comp1(:,1)-comp1(:,2))./abs(comp1(:,2));
%% Tercer Punto
% Orden n=5
[Hn5,b5]=hilbert(5); %Se genera una matriz 'Hn5' de Hilbert de dimensi�n n=5 y un vector columna 'b5' de  5x1
[H5,u5]=MatrizLU(Hn5);% Ejecuta el algoritmo de matriz LU. 'H5' matriz LU, 'u5' matriz de pivotes parciales calculados sobre 'Hn5'
[L5,U5,P5]=factorizacionLU(H5,u5);% Ejectuta el algoritmo para obtener las matrices L,U,P. 'L5' matriz diagonal inferior, 'U5' matriz diagonal superior, 'P5' matriz de permutaci�n.
y5=L5\(P5*b5);%Se soluciona el sistema Ly=Pb. 'L5'\('P5'*'b5')=inv(L)*P*b. 'y5' es el vector soluci�n.
x5=U5\y5; %Se soluciona el sistema Ux=y. 'U5'\'y5'=inv(U)*y. 'x5' es el vector soluci�n.
xt5=invhilb(5)*b5; %Soluci�n t�orica {Hn}x=b. 'Hnt5'\'b5'=inv(Hn)*b. 'xt5' es la soluci�n te�rica.
comp2=[x5 xt5]; %Se unen en una matriz los vectores columnas con la soluc�n mediante el algoritmo y la te�rica respectivamente
error2=abs(comp2(:,1)-comp2(:,2))./abs(comp2(:,2));

% Orden n=15
[Hn15,b15]=hilbert(15); %Se genera una matriz 'Hn15' de Hilbert de dimensi�n n=15 y un vector columna 'b15' de  15x1
[H15,u15]=MatrizLU(Hn15);% Ejecuta el algoritmo de matriz LU. 'H15' matriz LU, 'u15' matriz de pivotes parciales calculados sobre 'Hn15'
[L15,U15,P15]=factorizacionLU(H15,u15);% Ejectuta el algoritmo para obtener las matrices L,U,P. 'L15' matriz diagonal inferior, 'U15' matriz diagonal superior, 'P15' matriz de permutaci�n.
y15=L15\(P15*b15);%Se soluciona el sistema Ly=Pb. 'L15'\('P15'*'b15')=inv(L)*P*b. 'y15' es el vector soluci�n.
x15=U15\y15;%Se soluciona el sistema Ux=y. 'U15'\'y15'=inv(U)*y. 'x15' es el vector soluci�n.
xt15=invhilb(15)*b15; %Soluci�n t�orica {Hn}x=b. 'Hnt15'\'b15'=inv(Hn)*b. 'xt15' es la soluci�n te�rica.
comp3=[x15 xt15]; %Se unen en una matriz los vectores columnas con la soluc�n mediante el algoritmo y la te�rica respectivamente
error3=abs(comp3(:,1)-comp3(:,2))./abs(comp3(:,2));
% Orden n=53

[Hn53,b53]=hilbert(53);%Se genera una matriz 'Hn53' de Hilbert de dimensi�n n=53 y un vector columna 'b53' de  53x1
[H53,u53]=MatrizLU(Hn53);% Ejecuta el algoritmo de matriz LU. 'H53' matriz LU, 'u53' matriz de pivotes parciales calculados sobre 'Hn53'
[L53,U53,P53]=factorizacionLU(H53,u53);% Ejectuta el algoritmo para obtener las matrices L,U,P. 'L53' matriz diagonal inferior, 'U53' matriz diagonal superior, 'P53' matriz de permutaci�n.
y53=L53\(P53*b53);%Se soluciona el sistema Ly=Pb. 'L53'\('P53'*'b53')=inv(L)*P*b. 'y53' es el vector soluci�n.
x53=U53\y53;%Se soluciona el sistema Ux=y. 'U53'\'y53'=inv(U)*y. 'x53' es el vector soluci�n.
xt53=invhilb(53)*b53; %Soluci�n t�orica {Hn}x=b. 'Hnt53'\'b53'=inv(Hn)*b. 'xt53' es la soluci�n te�rica.
comp4=[x53 xt53]; %Se unen en una matriz los vectores columnas con la soluc�n mediante el algoritmo y la te�rica respectivamente
error4=abs(comp4(:,1)-comp4(:,2))./abs(comp4(:,2));
%%
function [L,U,P]=factorizacionLU(A,u1) %Genera las matrices L,U,P a partir de una matriz 'A' de la forma LU y una matriz de pivotes parciales 'u1'
%'A' es una matriz cuadrada forma LU INPUT
%'u1' es el vector de pivotes parciales INPUT
%'L' es una matriz cuadrada diagonal inferior OUTPUT
%'U' es una matriz cuadrada diagonal superior OUTPUT
%'P' es una matriz cuadrada de permutaci�n OUTPUT.
n=size(A,1); %se obtiene la dimensi�n de 'A'; 
L=ones(n); L=tril(L); %L se inicializa como una matriz diagonal de dimensi�n 'n' y unos en sus componentes
U=[];%U se inicializa como una matriz vacia
for i=1:n
    U(i,i:n)=A(i,i:n);% 'U' se le asigna los elementes de la matriz 'A' en la fila i desde la columna i+1 hasta 'n'.
    L(i+1:n,i)=A(i+1:n,i); %'L' se le asigna los elementos de la matriz 'A' en desde la fila i+1 hasta 'n' en la columna i.
end
P=eye(n);
for k=1:1:length(u1)
    P([k u1(k)],: ) = P([u1(k) k],:); %Se intercambia la fila k por la fila u(k)
end
end

%%
%Matriz de Hilbert
function [H,b] = hilbert(n) %Genera la matriz de Hilbert 'H' para un orden 'n' ingresado por par�metro y una matriz de resultados 'b'
%'H' es una matriz cuadrada de dimensi�n 'n' OUTPUT
%'b' es un vector columna de dimensi�n 1xn OUTPUT
%'n' es el orden o dimensi�n de las matrices INPUT
H=zeros(n); %Se inicializa 'H' como una matriz de ceror de dimensi�n nxn.
b=zeros(n,1);%Se inicializa 'b' como un vector columna de nx1
for i=1:n
    for j=1:n
       H(i,j)=1/(i+j-1); %Se modifica la posici�n (i,j) de 'H' por la expresi�n que depende de su posici�n.
       b(i)=b(i)+H(i,j);%Se modicia el i-�simo elemento del vector 'b' por la suma del mismo elemento en la posici�n i m�s el elemento en la posici�n (i,j) de la matriz 'H'
    end
end  
end

%%
function [A,u] = MatrizLU(A) %Genera las matrices de la forma LU y la matriz de pivotes parciales para una matriz 'A' dada.
%'A' es una matriz cuadrada e invertible INPUT
%'A' es la modificaci�n sobre ella misma que retorna en la forma LU OUTPUT
n=size(A,1); %se obtiene la dimensi�n de 'A'
for k=1:n-1
    u(k)=find(abs(A(k:n,k))==max(max(abs(A(k:n,k)))))+n-length(k:n); %Busca en la columna k-�sima desde la fila k hasta 'n' el indice del elemento m�ximo de la matriz 'A'.
 
    A([k u(k)],: ) = A([u(k) k],:); %Se intercambia la fila k por la fila u(k)
    if A(k,k)~=0
        p=k+1:n;
        A(p,k)=A(p,k)/A(k,k); %Ejecuta la operaci�n sobre la sub-columna k.
        A(p,p)=A(p,p)-A(p,k)*A(k,p); %Ejecuta la operaci�n sobre la sub-matriz A(p:n,p:n)
    else 
    end
end

end