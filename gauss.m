function U = gauss(A,F)
  N = size(A)(1);
  % On réalise les opérations sur les lignes
  for k = [N:-1:2]
    A(k,:) -= A(k-1,:);
    F(k) -= F(k-1);
  endfor
  
  % On réalise les opérations sur les colonnes
  lams = zeros(1,length([N:-1:2]));
  for k = [N:-1:2]
    lam = A(k,k-1)/A(k,k);
    lams(k-1) = lam;
    A(:,k-1) -= lam .* A(:,k);
  endfor
  %display(["A = "; num2str(A)]);
  U = solve_gauss(A,F);
  U = transpose(U);
endfunction
