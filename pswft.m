function I=pswft(f,N)
  
  % f : fonction dont on veut calculer le produit scalaire avec les polyomes de Tchebychev
  % N : nombre de produits scalaires calculés
  %
  % I : vecteur contenant les valeurs successives des produits scalaires
 
  tt=linspace(0,2*pi,N);
  w=f(tt);
  fft=fft(w);
  m=length(tt);
  I=pi*real(fft)./N;
endfunction
