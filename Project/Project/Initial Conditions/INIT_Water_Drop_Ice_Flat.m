function L = INIT_Water_Drop_Ice_Flat(x,P)
h0=1.6e-9; Af= 17/P.k1; xwf = 450/P.k1;
Nx=length(x); xL=x(end);
% Set the initial condition as a (symmetric) Gaussian function
Llv_init= P.dB + h0 + Af*ones(Nx,1).*(exp(-((x-xL/2)/xwf).^2))';
%The ice is just a flat sheet
Lsl_init=P.dB*ones(Nx,1);
L=[Llv_init;Lsl_init];
end

