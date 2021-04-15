function U = reorder(V)
  l = length(V);
  U = zeros(l,1);
  K = floor(l/2);
  for k = [1:K]
    U(2*k-1) = V(k);
    U(2*k) = V(K+k);
  endfor
endfunction
