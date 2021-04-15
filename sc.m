function I=sc(f,n,N)
  t=linspace(0,2*pi,N);
  g=f(t).*cos(n*t);
  I=pi*sum(g)/N;
endfunction
