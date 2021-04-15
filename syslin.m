function [M1,M2,G1,G2] = syslin(alpha,gamma,f,N)
  % On commence par créer le système linéaire pour lequel on n'a pas pris en compte la parité
  A = zeros(N,N);
  B = zeros(N,N);
  F = zeros(N,1);
  for l = [0:N-1]
    for k = [0:N-1]
      if mod(l,2) == mod(k,2)
        if mod(l,2) == 0
          if l==k
            A(l+1,k+1) += 3*pi/3;
          else 
            A(l+1,k+1) += pi;
          endif
          if l < k
            B(l+1,k+1) += -pi/2 * l**2 * k;
          else
            B(l+1,k+1) += -pi/2 * k**3;
          endif
        else
          if l==k
            A(l+1,k+1) += pi;
          else 
            A(l+1,k+1) += pi/2;
          endif
          if l < k
            B(l+1,k+1) += -pi/2 * (l**2 - 1) * k;
          else
            B(l+1,k+1) += -pi/2 * (k**2 -1)*k;
          endif
        endif
      endif
    endfor
    F(l+1) = pswft(f,l,80);
  endfor
  A = -alpha* B + gamma * A;
  %display(["A = "; num2str(A)]);
  % On modifie à présent le système pour prendre en compte la parité
  K = floor(N/2);
  M1 = zeros(K,K);
  M2 = zeros(K,K);
  G1 = zeros(K,1);
  G2 = zeros(K,1);
  for l = [1:K]
    for k = [1:K]
      M1(l,k) = A(2*l-1,2*k-1);
      M2(l,k) = A(2*l,2*k);
    endfor
    G1(l) = F(2*l-1);
    G2(l) = F(2*l);
  endfor
endfunction
