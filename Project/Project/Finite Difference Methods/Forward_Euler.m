function L = Forward_Euler(odefun,dt,Nt,L,P,D1,D2)
for i = 1:Nt-1
    L = L + dt*odefun(dt,L,P,D1,D2);
end
end

