clear; clc; close all;
p = 1;
f = @(x)(cos(pi*x));
alpha = 1;
gamma = 10;
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

U = [U1;U2];
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
  Tk = @(t)(cos((k+1)*t));
  Tj = @(t)(cos(j*t));
  norm_k = sqrt(scwft(Tk,k+1,J) + scwft(Tj,j,J));
  Phi(k,:) /= norm_k;
  if k < 6
    plot(x,Phi(k,:),"Displayname",num2str(k-1))
  endif
endfor
legend
hold off

sol = transpose(U) * Phi;

[M1,M2,G1,G2] = syslin(alpha,gamma,f,N);
U1 = gauss(M1,G1);
U2 = gauss(M2,G2);
U = [U1;U2];
U = reorder(U);
sol = transpose(U) * Phi;


[A,F] = syslin_df(gamma,f,0,0,J+1);
sol_df = A\F;

figure
hold on
plot(x,sol,"Displayname","Solution")
plot(x,sol_df,"Displayname","Solution DF")
legend
hold off

