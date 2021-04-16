function u=sedoci(x,alpha,dalpha,beta,gamma,f,ci)
% Calcule une approximation de la solution de l'équation différentielle
%   -(alpha u'(x))' + beta(x) u'(x) + gamma(x) u(x) = f(x)
% sur l'intervalle [a,b] 
% sous les conditions initiales
%   u(a) = ci(1)
%   u'(a) = ci(2)
% par la méthode d'Euler
%
% Paramètres d'entrée :
% x : vecteur ligne contenant les noeuds de discretisation de l'intervalle [a,b] 
% alpha, beta, gamma : fonctions coefficients de l'edo
% f : second membre de l'edo
% dalpha : dérivée de la fonction alpha
% ci : vecteur colonne contenant les données initiales [u(a);u'(a)]
%
% Paramètres de sorties :
% u : la solution calculée aux noeuds de discrétisation 
%
% Auteur : Stéphane Balac - UFR de Mathématiques - Université de Rennes 1
% Décembre 2015

A=@(x) [0,1;gamma(x)./alpha(x),(beta(x)-dalpha(x))./alpha(x)];
B=@(x) [0; -f(x)./alpha(x)];

if (size(ci,1)==1), ci=transpose(ci);end; % on s'assure que le vecteur est bien en colonne
Y=zeros(2,length(x)); % déclaration du tableau
Y(:,1)=ci; % initialisation

for i=1:length(x)-1
	h=x(i+1)-x(i); % pas courant
	Y(:,i+1)=Y(:,i)+h*(A(x(i))*Y(:,i)+B(x(i))); % schéma d'Euler
end
u=Y(1,:);
end