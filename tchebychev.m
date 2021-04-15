% Retourne les N premiers polynomes de tchebychev calculés sur x
function T = tchebychev(N,x)
  l = length(x);
  T = zeros(N,l);
  T(1,:) = ones(1,l);
  T(2,:) = x;
  for k = [3:N]
    T(k,:) = 2*x .* T(k-1,:) - T(k-2,:); 
  endfor
endfunction
