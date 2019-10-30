clc
clear
%%

m=1000;
n=100;
p=m;
A=1000*rand(m);
B=100*rand(p,n);

%Orden i,j,k
C=zeros(m,n);
for i=1:m
    for j=1:n
        for k=1:p
            C(i,j)=C(i,j)+A(i,k)*B(k,j);
        end
    end
end
t3=toc

%Orden i,k,j

C1=zeros(m,n);
for i=1:m
    for k=1:p
        for j=1:n
            C1(i,j)=C1(i,j)+A(i,k)*B(k,j);
        end
    end
end

%Orden k,j,i

C2=zeros(m,n);
for k=1:p
    for j=1:n
        for i=1:m
            C2(i,j)=C2(i,j)+A(i,k)*B(k,j);
        end
    end
end

%Orden k,i,j

C3=zeros(m,n);
for k=1:p
    for i=1:m
        for j=1:n
            C3(i,j)=C3(i,j)+A(i,k)*B(k,j);
        end
    end
end

%Orden j,i,k

C4=zeros(m,n);
for j=1:n
    for i=1:m
        for k=1:p
            C4(i,j)=C4(i,j)+A(i,k)*B(k,j);
        end
    end
end

%Orden j,k,i

C5=zeros(m,n);
for j=1:n
    for k=1:p
        for i=1:m
            C5(i,j)=C5(i,j)+A(i,k)*B(k,j);
        end
    end
end


%%
% 2 Recorridos
%Permutacion (i,j,k)
tic
C2ij=zeros(m,n);
for i=1:m
    for j=1:n
        C2ij(i,j)=C2ij(i,j)+A(i,:)*B(:,j);
    end
end
toc
t2=toc;
%Permutacion (i,k,j)
C2ik=zeros(m,n);
for i=1:m
    for k=1:n
        C2ik(i,k)=C2ik(i,k)+A(i,:)*B(:,k);
    end
end
%Permutacion (k,j,i)
C2kj=zeros(m,n);
for k=1:m
    for j=1:n
        C2kj(k,j)=C2kj(k,j)+A(k,:)*B(:,j);
    end
end

%Permutacion (k,i,j)
C2ki=zeros(m,n);
for k=1:m
    for i=1:n
        C2ki(k,i)=C2ki(k,i)+A(k,:)*B(:,i);
    end
end

%Permutacion (j,k,i)
C2jk=zeros(m,n);
for j=1:m
    for k=1:n
        C2jk(j,k)=C2jk(j,k)+A(j,:)*B(:,k);
    end
end

%Permutacion (j,i,k)
C2ji=zeros(m,n);
for i=1:m
    for j=1:n
        C2ji(i,j)=C2ji(i,j)+A(i,:)*B(:,j);
    end
end

%%
% 1 Recorrido
%Permutacion (i,j,k),(i,k,j)
tic
C1i=zeros(m,n);
for i=1:m
    C1i(:,:)=C1i(:,:)+A(:,i)*B(i,:);
end
toc
%Permutacion (k,j,i),(k,i,j)
C1k=zeros(m,n);
for k=1:m
    C1k(:,:)=C1k(:,:)+A(:,k)*B(k,:);
end
%Permutacion (j,k,i),(j,i,k)
C1j=zeros(m,n);
for j=1:m
    C1j(:,:)=C1j(:,:)+A(:,j)*B(j,:);
end
t1=toc;

