function w = spectral_D1(v,P)
Nx=length(v);
v_hat = fft(v);
%the wave frequencies
k = [0:Nx/2-1, 0, -Nx/2+1:-1]';
w_hat = 1i*k.*v_hat;
scale_factor=2*pi*P.k1/2500;
w=real(ifft(w_hat))*scale_factor;
end

