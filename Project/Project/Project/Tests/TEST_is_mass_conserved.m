function mass = TEST_is_mass_conserved(tolerance, time_method,odefun,t,dt,x,L,P,D1,D2)
        mass_before = sum(L(1,:))*P.rhol + sum(L(2,:))*P.rhos;
        L=time_method(odefun,dt,round(t*P.tau/dt),L,P,D1,D2);
        mass_after = sum(L(1,:))*P.rhol + sum(L(2,:))*P.rhos;
        mass=[mass_before,mass_after];
        relative_mass = (mass_before-mass_after)/mass_before;
        if(abs(relative_mass))<tolerance
            is_conserved = true;
        else 
            is_conserved = false;
        end
end
