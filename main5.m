%%% Script pour l'exercice 5 - transformee de Fourier rapide

clear; clc; close all;

p=1; n=6; N=255;
fcos = @(t) ((t.^p).*((2*pi.-t)).^p);

%%% Trace de fcos

figure();
fplot(fcos, [0, 2*pi], "DisplayName", "fcos");
xlabel("tet");
ylabel("f(cos(tet))");
title(["Trace de fcos pour p=",num2str(p)]);
J = 20;
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
plot(N,log(errs));
title("Erreur relative");
xlabel("N");
ylabel("Erreur");

%%% Temps de calcul compare a une methode d'integration naive

p=4; N=255;
tft=[];
tint=[];
Iint=[];
Ifft=[];
for n=0:N
  fcos = @(tet) (tet.^p).*((2*pi.-tet).^p);
  ts=cputime;
  Ifft(n+1)=scwft(fcos,n,N);
  te=cputime;
  tft(n+1)=te-ts;
  fun=@(x) fcos(x).*cos(n.*x);
  ts=cputime;
  Iint(n+1)=integral(fun,0,2*pi)/2;
  te=cputime;
  tint(n+1)=te-ts;
endfor

disp(["temps total de calcul par transformée de fourier : ",num2str(sum(tft))," s"]);
disp(["temps total de calcul par integration numerique : ",num2str(sum(tint))," s"]);