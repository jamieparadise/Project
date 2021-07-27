function L = INIT_Water_Flat_Ice_Square(x,P)
h0=2.5*P.dB;%1.6e-9; 
Ai= 1; 
Nx=length(x);xL= x(end); xwi=xL/16;
%The water is just a flat 
Llv_init=(P.dB+h0)*ones(Nx,1);
Lsl_init=P.dB+(Ai/2)*P.dB*(tanh(P.k1*(x-(xL-xwi)/2)/10)-tanh(P.k1*(x-(xL+xwi)/2)/10))';
L=[Llv_init;Lsl_init];
end

