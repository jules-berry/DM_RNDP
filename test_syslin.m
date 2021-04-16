% Script qui permet de tester la méthode
clear; clc; close all;
p = 1;
f = @(x)(x);
alpha = 1;
gamma = 1;
N = 20;
[M1,M2,G1,G2] = syslin(alpha,gamma,f,N);
%display(["M1 = "; num2str(M1)])
%display(["M2 = "; num2str(M2)])
%display(["G1 = "; num2str(G1)])
%display(["G2 = "; num2str(G2)])

U1 = gauss(M1,G1);
U2 = gauss(M2,G2);

%display(["U1 = "; num2str(U1)])
%display(["U2 = "; num2str(U2)])

V1 = M1\G1;
V2 = M2\G2;

%display(["V1 = "; num2str(V1)])
%display(["V2 = "; num2str(V2)])

% On trouve que les deux méthodes de resolution donnent le meme resultat.

U = [V1;V2];
%display(["U = "; num2str(U)])

% On réordonne les éléments de U
U = reorder(U);
%display(["U = "; num2str(U)])

J = 250;
x = linspace(-1,1,J);
T = tchebychev(N+2,x);
figure
hold on
for k = [1:5]
  plot(x,T(k,:),"Displayname",num2str(k-1),"LineWidth",1.3)
endfor
title("Polynomes de Tchebychev")
xlabel("x")
ylabel("T_k(x)")
legend
hold off

Phi = zeros(N,J);
figure
hold on
for k = [1:N]
  j = mod(k-1,2);
  Phi(k,:) = T(k+2,:) - T(j+1,:);
  Tk = @(t)(cos((k+1)*t));
  Tj = @(t)(cos(j*t));
  norm_k = (abs(pswft(Tk,max(80,N+3))(k+3)) + abs(pswft(Tj,max(80,N+3))(j+1)));
  Phi(k,:) /= norm_k;
  if k < 6
    plot(x,Phi(k,:),"Displayname",num2str(k-1),"LineWidth",1.3)
  endif
endfor
xlabel("x")
ylabel("phi_k(x)")
title("Fonctions phi")
legend
hold off

[M1,M2,G1,G2] = syslin(alpha,gamma,f,N);
U1 = gauss(M1,G1);
U2 = gauss(M2,G2);
V1 = M1\G1;
V2 = M2\G2;
U = [V1;V2];
U = reorder(U);
sol = transpose(U) * Phi;


%sol_ex = @(x)(f(x)/gamma *(1 - cosh(sqrt(gamma) * x).* 1/cosh(sqrt(gamma))));
u1 = @(x)(2*exp(-x -1) - exp(x + 1) + x);
u2 = @(x)(1/2 * (exp(x + 1) - exp(-x-1)));
kappa = - u1(1)/u2(1) ;
sol_ex = @(x)(u1(x) + kappa * u2(x));

u1_exp = sedoci(x,@(xt)(alpha),@(t)(0),@(t)(0),@(t)(gamma),f,[0;0]);
u2_exp = sedoci(x,@(xt)(alpha),@(t)(0),@(t)(0),@(t)(gamma),f,[0;1]);
kappa = - u1_exp(J)/u2_exp(J);

sol_exp = u1_exp + kappa* u2_exp;

figure
hold on
plot(x,sol,"Displayname","Solution","LineWidth",1.5)
plot(x,sol_exp,"Displayname","Solution Tir","Color","black","LineWidth",1.5)
plot(x,sol_ex(x),"Displayname","Solution exacte","Color","red","LineWidth",1.5);
xlabel("x")
ylabel("u(x)")
title("Résultats pour -u'' + u = x ; u(-1) = u(1) = 0")
legend
hold off

