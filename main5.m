%%% Script pour l'exercice 5 - transformee de Fourier rapide

clear; clc; close all;

p=1; n=6; N=255;
fcos = @(t) ((t.^p).*((2*pi.-t)).^p);
%trac� de fcos
figure();
fplot(fcos, [0, 2*pi]);
xlabel("tet");
ylabel("f(cos(tet))");
title(["Trace de fcos pour p=",num2str(p)]);
J = 20
res=zeros(1,J);
K = 1:J;
coefs_reels = -2*pi./(K.^2);
coefs_reels = [2*pi**3/3,coefs_reels];
%resint=[];
for i=0:J
  res(i+1)=scwft(fcos,i,N);
 % resint(i)=sc(f,i,N);
endfor
figure();
hold on;
plot(0:J,res,"DisplayName","Coefs calculés");
plot(0:J,coefs_reels,"DisplayName","Coefs réels");
xlabel("n");
ylabel("pswft(fcos,n)");
title("Calcul des produits scalaires");
legend;
hold off;

%%% Convergence de la FFT

N = [20:5:200];
coefs = zeros(1,length(N));
errs = zeros(1,length(N));
for i=[1:length(N)]
  coefs(i) = scwft(fcos,5,N(i));
  errs(i) = abs(-2*pi/25-coefs(i))/(2*pi/25);
endfor

figure;
plot(N,errs);
title("Erreur relative");
xlabel("N");
ylabel("Erreur");

%%% Temps de calcul compare a une methode d'integration naive

p=4; n=5; N=255;
fcos = @(tet) (tet.^p).*((2*pi.-tet).^p);
g=@(x) fcos(x).*cos(n*x);