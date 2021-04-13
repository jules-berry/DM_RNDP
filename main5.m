%%% Script pour l'exercice 5 - transformee de Fourier rapide

clear; clc; close all;

p=1; n=6; N=255;
fcos = @(tet) (tet.^p).*((2*pi.-tet).^p);

%tracé de fcos
figure();
fplot(fcos, [-2*pi, 2*pi]);
xlabel("tet");
ylabel("f(cos(tet))");
title(["Trace de fcos pour p=",num2str(p)]);

res=[];
resint=[];
for i=1:N
  res(i)=scwft(fcos,n,i);
  resint(i)=sc(fcos,n,i);
endfor
figure();
plot(1:N,res);
xlabel("n");
ylabel("pswft(fcos,n)");
title("Calcul des produits scalaires");

%%% Convergence de la FFT




%%% Temps de calcul compare a une methode d'integration naive

p=4; n=5; N=255;
fcos = @(tet) (tet.^p).*((2*pi.-tet).^p);
g=@(x) fcos(x).*cos(n*x);