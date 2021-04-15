function [A,F]=syslin_df(gamma,f,u0,u1,J)
% Construction du système linéaire issu de la discrétisation du problème
% aux limites de Dirichlet
% -u''(x) + gamma(x) u(x) = f(x) pour tout x dans ]0,1[
% u(0)=u0   u(1)=u1
% par un schéma aux différences finies basé sur une approximation de la
% dérivée seconde par une différence divisée centrée d'ordre 2.
% L'intervalle [0,1] est subdivisé en J sous-intervalles.

h=1/J; % pas de la subdivision
x=linspace(-1,1,J-1); % noeuds de la subdivision
B=[-ones(1,J-1);2*ones(1,J-1);-ones(1,J-1)];
B=transpose(B)./h^2;
B(:,2)=B(:,2)+ gamma;
A=spdiags(B,[-1,0,1],J-1,J-1);
F=transpose(f(x));
F(1)=f(x(1))+u0/h^2;
F(J-1)=f(x(J-1))+u1/h^2;