function P = Load_Parameters(pv,T)
P.pv=pv; %pressure in Pa
P.T=T; %temp in K
P.Tc=P.T - 273.15; %temp in C

%Parameter values  from Sup table 1
P.C1 = 3.143e-3;%;
P.C2 =4.116e-2; %
P.k1 = 2.043e9;
P.k2 =2*P.k1;
P.q0 = 7.148e9;
P.alpha = 1.114; %should be 5.114
P.B = 7.785e-31;
P.f= 1.106;
P.a = 3.03e7;
P.b= 5.0e8;
%Parameter values  from Sup table 2
P.Tt = 273.16;
%P.Ttc = 0.1;
P.pt = 611.65;
% rhol = 55498; 
% rhos = 50888;
% rhov  = 0.2694;
P.delta_Hsv= 51059;
P.delta_Hlv = 45051;
P.delta_Hsl = 6008;
%Parameter values  from Sup table 3
P.R=8.31446261815324;
P.dB = 0.37e-9;%
%P.Rc= 461.52;
P.rhol = 55502 + 3.4549*P.Tc - 0.44461*P.Tc^2 + 0.0028885*P.Tc^3 -0.00031898*P.Tc^4;%
P.rhos = 50885 - 9.71*P.Tc - 0.03*P.Tc^2;%
P.plv=453;
P.rholv = P.plv/(P.R*P.T);
P.rhov = P.pv/(P.R*P.T);

P.gamma_sl=(28 + 0.25*P.Tc)*10^-3; 
P.gamma_lv=(75.7 - 0.1775*P.Tc)*10^-3; 
P.eta =(1.39e-4)*(P.T/225-1)^-1.64; 
P.u = ((P.dB/(2*pi))*P.rhos*P.delta_Hsl*2)/P.Tt;%1.3e-4
P.klv =(3.4e-10)*P.rholv*P.T^(-1/2); %should be *10^-3
P.ksl=6.4*P.klv;% 
%Parameter values  from Sup Note 6

P.delta_psl=((P.rhos-P.rhol)*P.R*P.T*log(P.pv/P.pt)+(P.rhol*P.delta_Hlv-P.rhos*P.delta_Hsv)*(P.T-P.Tt)/P.Tt);%
P.delta_plv=(P.rhol*P.R*P.T*log(P.pv/P.pt)-P.rhol*P.delta_Hlv*(P.T-P.Tt)/P.Tt);%
P.delta_pk=(P.rhos*P.ksl*P.delta_psl - P.rhol*P.klv*P.delta_plv)/(P.rhos*P.ksl+P.rhol*P.klv); %this is wrong
P.tau = 3.8136e-11;% (3*P.eta)/(P.k1*P.gamma_lv)=3.8136e-11;%  using this formula. supposed to be 0.11e-9 in the paper
% More parameters
P.delta_rho = P.rhos-P.rhol;
P.qz=2*pi/P.dB;
end

