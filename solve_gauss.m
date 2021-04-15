function U = solve_gauss(A,F)
  N = size(A)(1);
  U = zeros(1,N);
  U(N) = F(N)/A(N,N);
  for k = [N-1:-1:1]
    U(k) = (F(k) - sum(A(k,k+1:end).*U(k+1:end)))/A(k,k);
  endfor
endfunction
