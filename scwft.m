function I=scwft(f,n,N)
  tt=linspace(0,2*pi,N);
  g=f(tt).*cos(n*tt);
  
  fourier=fft(g);
  m=length(tt);
  fourier=fourier./m;
  %display(num2str(m))
  I=2*pi*fourier(1);
  J=1:m;
  
  for k=1:floor(m/2) %attention peut etre +1
    S = sum(cos(2*pi*k.*J./m));
    I = I + (2*pi)/m * (fourier(1+k) + fourier(m - k)) * S;
  endfor
  I = real(I/2);
endfunction
