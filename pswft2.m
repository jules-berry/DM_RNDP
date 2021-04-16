function I = pswft2(f,N)
  % f : fonction dont on veut calculer le produit scalaire avec les polyomes
  % de Tchebychev
  % N : nombre de produits scalaires calculés
  %
  % I : vecteur contenant les valeurs successives des produits scalaires
  fcos = @(x)(f(cos(x)));
  I = pswft(fcos,N);
endfunction
