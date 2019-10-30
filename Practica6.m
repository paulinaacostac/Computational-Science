%%
clc
clear
%%
A1=[4 1 -1; -1 3 1; 2 2 5]; b1=[5;-4;1];
A2=[-2 1 1/2;1 -2 -1/2;0 1 2]; b2=[4;-4;0];
A3=[4 1 -1 1;1 4 -1 -1;-1 -1 5 1;1 -1 1 3]; b3=[-2;-1;0;1];
w=1:0.01:1.5;

[xp,kp,difp]=GaussSeidel(A3,b3,1,100,1e-6);

iter1=zeros(length(w),3); %Vector resultado variando parametros de w, primera columna w, segunda columna N, tercera columna dif
iter2 = zeros(length(w),3);%Vector resultado variando parametros de w, primera columna w, segunda columna N, tercera columna dif
iter3 = zeros(length(w),3); %Vector resultado variando parametros de w, primera columna w, segunda columna N, tercera columna dif

for i=1:length(w)
    
    [x,k,dif]=GaussSeidel(A1,b1,w(i),100,1e-6);
    iter1(i,1)=w(i);
    iter1(i,2)=k;
    iter1(i,3)=dif;
    [x1,k1,dif1]=GaussSeidel(A2,b2,w(i),100,1e-6);
    iter2(i,1)=w(i);
    iter2(i,2)=k1;
    iter2(i,3)=dif1;
    [x2,k2,dif2]=GaussSeidel(A3,b3,w(i),100,1e-6);
    iter3(i,1)=w(i);
    iter3(i,2)=k2;
    iter3(i,3)=dif2;
 
end
T1 = array2table(iter1,...
    'VariableNames',{'omega','N','dif'});

T2 = array2table(iter2,...
    'VariableNames',{'omega','N','dif'});

T3 = array2table(iter3,...
    'VariableNames',{'omega','N','dif'});


figure
plot(iter1(:,1),iter1(:,2))
title('Numero de iteraciones en función de \omega sistema A')
xlabel('\omega')
ylabel('Número de iteraciones N')
set(gcf,'color','w')


figure
plot(iter2(:,1),iter2(:,2))
title('Numero de iteraciones en función de \omega sistema B')
xlabel('\omega')
ylabel('Número de iteraciones N')
set(gcf,'color','w')

figure
plot(iter3(:,1),iter3(:,2))
title('Numero de iteraciones en función de \omega sistema C')
xlabel('\omega')
ylabel('Número de iteraciones N')
set(gcf,'color','w')

n=80;
A=matrizA(n);
b=pi*ones(n,1);

wp=[1 1.5 1.9]';
iter4=zeros(length(wp),3);

for i=1:length(wp)
    [x1,k1,dif1]=GaussSeidel(A,b,wp(i),200,1e-6);
    iter4(i,1)=wp(i);
    iter4(i,2)=k1;
    iter4(i,3)=dif1;
end

T4 = array2table(iter4,...
    'VariableNames',{'omega','N','dif'});


%%
function [A]=matrizA(n)
A=zeros(n);
for i=1:n
    for j=1:n
        if i==j
            A(i,j)=4;
        elseif (j==i-1) ||(j==i+1)
            A(i,j)=-1;   
        end
    end
end


end

%%
function [x,k,dif]=GaussSeidel(A,b,w,N,eps)

k=1;
dif=1;
n=length(A);
x0=zeros(n,1);
x=zeros(n,1);
while k<=N && dif>eps
    
        for i=1:n
            x(i)=(1-w)*x0(i)+(w/A(i,i))*(-A(i,1:i-1)*x(1:i-1)- A(i,i+1:n)*x0(i+1:n)+b(i));
        end
  
    dif=max((abs(x0-x)));
    k=k+1;
    x0=x;
    
end

end