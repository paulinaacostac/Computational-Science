%%
%P. Acosta, J. Timana. 2019-1 Universidad de los Andes: Curso de
%computación científica
clc
clear
%%
%Primer punto
format shortE
%%
%Segundo punto
r=1/3; %argumento serie geometrica
n=10; %numero de iteraciones
a=1; %constante

%%%Cálculo mediante formula (exacto)%%%
SerieGeo1=a*(1-r^n)/(1-r);

%%%Cálculo mediante método aproximado%%%

SerieGeo2=0;
for i=1:n
SerieGeo2=SerieGeo2+a*r^(i-1);
end

%%%Comparación%%%
Comp= SerieGeo1==SerieGeo2;
%%
cos1=cos(2*pi);
cos2=cos(pi/2);
numMax = realmax;
numMin = realmin;
