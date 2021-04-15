clear; clc; close all;
p = 1;
fcos = @(t) ((t.^p).*((2*pi.-t)).^p);
alpha = 1;
gamma = 1;
N = 10;
[M1,M2,G1,G2] = syslin(alpha,gamma,fcos,N);
display(["M1 = "; num2str(M1)])
display(["M2 = "; num2str(M2)])
display(["G1 = "; num2str(G1)])
display(["G2 = "; num2str(G2)])

U1 = gauss(M1,G1);
U2 = gauss(M2,G2);

display(["U1 = "; num2str(U1)])
display(["U2 = "; num2str(U2)])

V1 = M1\G1;
V2 = M2\G2;

display(["V1 = "; num2str(V1)])
display(["V2 = "; num2str(V2)])

% On trouve que les deux méthodes de résolution donnent le meme résultat.

U = [U1;U2];
display(["U = "; num2str(U)])

% On réordonne les éléments de U
U = reorder(U);
display(["U = "; num2str(U)])

J = 250;
x = linspace(0,1,J);
T = tchebychev(N+2,x);
Phi = zeros(N,J);
for k = [1:N]
  j = mod(3,2);
  Phi(k,:) = T(k+2,:) - T(j,:);
endfor

sol = sum(U .* Phi);

figure
hold on
plot(x,sol,"Displayname","Solution")
legend
hold off

