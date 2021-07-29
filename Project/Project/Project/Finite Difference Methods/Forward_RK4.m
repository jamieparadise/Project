function L = Forward_RK4(odefun,dt,Nt,L,P,D1,D2)
for i = 1:Nt-1
    K1 = odefun((i-1)*dt,L,P,D1,D2);
    K2 = odefun((i-.5)*dt,L+dt*K1/2,P,D1,D2);
    K3 = odefun((i-.5)*dt,L+dt*K2/2,P,D1,D2);
    K4 = odefun(i*dt,L+dt*K3,P,D1,D2);
    L = L + dt*(K1/6+K2/3+K3/3+K4/6);
end
end