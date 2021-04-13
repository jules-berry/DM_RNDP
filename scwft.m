function I=scwft(f,n,N)
  tt=linspace(0,2*pi,N);
  g=f(tt).*cos(n*tt);
  
  fft=fft(g);
  m=length(tt);
  fft=fft./m;
  
  I=2*pi*fft(1);
  J=1:m;
  
  for k=1:floor(m/2) %attention peut etre +1
    S = sum(cos(2*pi*k.*J./m));
    I = I + (2*pi)/m * (fft(k) + fft(floor(m/2) + k)) * S;
  endfor
  I = real(I/2);
endfunction
