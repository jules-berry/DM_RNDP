clear; clc; close all;
p = 1;
f = @(x)(1);
alpha = 1;
gamma = 1;
N = 100;
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

% On trouve que les deux méthodes de résolution donnent le meme résultat.

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
  plot(x,T(k,:),"Displayname",num2str(k-1))
endfor
title("Polynomes de tchebychev")
legend
hold off

Phi = zeros(N,J);
figure
hold on
for k = [1:N]
  j = mod(k-1,2);
  Phi(k,:) = T(k+2,:) - T(j+1,:);
  Tk = @(t)(cos((k+1) * acos(t)));
  Tj = @(t)(cos(j*acos(t)));
  norm_k = (pswft(Tk,k+1,J) + pswft(Tj,j,J));
  Phi(k,:) /= norm_k;
  if k < 6
    plot(x,Phi(k,:),"Displayname",num2str(k-1))
  endif
endfor
legend
hold off

[M1,M2,G1,G2] = syslin(alpha,gamma,f,N);
U1 = gauss(M1,G1);
U2 = gauss(M2,G2);
V1 = M1\G1;
V2 = M2\G2;
U = [V1;V2];
U = reorder(U);
sol = zeros(1,length(x));
for k = [0:N-1]
  j = mod(k,2);
  Tk = @(t)(cos((k+2)*acos(t)));
  Tj = @(t)(cos(j*acos(t)));
  norm_k = (pswft2(Tk,k+2,J) + pswft2(Tj,j,J));
  sol += U(k+1) * (Tk(x) - Tj(x))./norm_k;
endfor


[A,F] = syslin_df(gamma,f,0,0,J+1);
sol_df = A\F;

sol_ex = @(x)(f(x)/gamma *(1 - cosh(sqrt(gamma) * x).* 1/cosh(sqrt(gamma))));

figure
hold on
plot(x,sol,"Displayname","Solution")
plot(x,sol_df,"Displayname","Solution DF")
plot(x,sol_ex(x),"Displayname","Solution exacte");
legend
hold off

