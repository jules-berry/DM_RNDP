function I = pswft2(f,n,N)
  J = [0:N-1];
  w = f(cos(2*pi*J./N));
  e = exp(-1i*pi*n*J./N);
  I = pi/N* real(sum(w.*e));
endfunction
