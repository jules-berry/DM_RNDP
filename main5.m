%%% Script pour l'exercice 5 - transformee de Fourier rapide %%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%%% Trace de fcos ==============================================================
p=1;
fcos = @(t) ((t.^p).*((2*pi.-t)).^p);
figure();
fplot(fcos, [0 , 2*pi], "DisplayName", "fcos");
xlabel("tet");
ylabel("f(cos(tet))");
title(["Trace de fcos pour p=",num2str(p)]);


%%% Trace des coefficients  ====================================================
J = 20;  N=85;
res=zeros(1,J);
K = 1:J;
coefs_reels = -2*pi./(K.^2);
coefs_reels = [2*pi**3/3,coefs_reels];
ps = pswft(fcos,N);
display(num2str(size(ps)));
%resint=[];
for i=0:J
  res(i+1)=ps(i+1);
endfor

figure();
hold on;
plot(0:J,res,"DisplayName","Coefs calcul�s");
plot(0:J,coefs_reels,"DisplayName","Coefs r�els");
xlabel("n");
ylabel("pswft(fcos,n)");
title("Calcul des produits scalaires");
legend;
hold off;


%%% Convergence de la FFT  =====================================================
p=1; n=5;
fcos = @(t) ((t.^p).*((2*pi.-t)).^p);
N = [20:5:100];
coefs = zeros(1,length(N));
errs = zeros(1,length(N));
for i=[1:length(N)]
  coefs(i) = pswft(fcos,N(i))(n+1);
  errs(i) = abs(-2*pi/(n^2)-coefs(i))/(2*pi/(n^2));
endfor

upto=length(N)-3;
P=polyfit(N(1:upto), log(errs)(1:upto),1);
disp(["coefficient directeur par interpolation : ",num2str(P(1))]);

figure;
hold on;
plot(N,[log(errs); P(2).+P(1).*N]);
legend(["log(err)"; "droite d'interpolation"]);
title(["Erreur relative pour (p,n)=(", num2str(p),",", num2str(n),")"]);
xlabel("N");
ylabel("log(Erreur)");
hold off;


%%% Trace de fcos ==============================================================
p=4;
fcos = @(t) ((t.^p).*((2*pi.-t)).^p);
figure();
hold on;
fplot(fcos, [0 , 2*pi], "DisplayName", "fcos");
xlabel("tet");
ylabel("f(cos(tet))");
title(["Trace de fcos pour p=",num2str(p)]);
hold off;

%{
%%% Convergence de la FFT  =====================================================
p=4; n=0;
fcos = @(t) ((t.^p).*((2*pi.-t)).^p);
N = [20:5:300];
coefs = zeros(1,length(N));
errs = zeros(1,length(N));

fun=@(x) fcos(x).*cos(n.*x);
valex=256*pi/(2*315); %pswft(fcos,85)(n+1); %a remplacer par la valeur exacte
for i=[1:length(N)]
  coefs(i) = pswft(fcos,N(i))(n+1);
  errs(i) = abs((valex-coefs(i))/(valex));
endfor

upto=length(N)-3;
P=polyfit(N(1:upto), log(errs)(1:upto),1);
disp(["coefficient directeur par interpolation : ",num2str(P(1))]);

figure;
plot(N,[log(errs); P(2).+P(1).*N]);
legend(["log(err)"; "droite d'interpolation"]);
title(["Erreur relative pour (p,n)=(", num2str(p),",", num2str(n),")"]);
xlabel("N");
ylabel("log(Erreur)");
%}

%%% Temps de calcul compare a une methode d'integration naive ==================
p=4; N=255;
Iint=zeros(1,N);
Ifft=zeros(1,N);

ts=cputime;
fcos = @(tet) (tet.^p).*((2*pi.-tet).^p);
Ifft=pswft(fcos,N);
te=cputime;
tft=te-ts;

ts=cputime;
for n=0:N
  fun=@(x) fcos(x).*cos(n.*x);
  Iint(n+1)=integral(fun,0,2*pi)/2;
endfor
te=cputime;
tint=te-ts;

disp(["temps total de calcul par transform�e de fourier : ",num2str(tft)," s"]);
disp(["temps total de calcul par integration numerique : ",num2str(tint)," s"]);
