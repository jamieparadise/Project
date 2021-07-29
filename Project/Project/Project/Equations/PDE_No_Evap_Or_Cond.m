function f = PDE_No_Evap_Or_Cond(t,L,P,D1,D2)
    Llv=L(:,1); Lsl=L(:,2);
    h=Llv-Lsl; % h is also temporary
    g_prime = -P.C2*P.k2*exp(-P.k2*h) + P.C1*exp(-P.k1*h).*(P.q0*sin(P.q0*h + P.alpha) + P.k1*cos(P.q0*h + P.alpha)) ...
    + 3*P.B*(h.^-4).*(1 - P.f*exp(-P.a*h)-(1-P.f)*exp(-P.b*h)) - P.B*(h.^-3).*( P.a*P.f*exp(-P.a*h)+P.b*(1-P.f)*exp(-P.b*h));
    DODLsl = P.u*P.qz*sin(P.qz*Lsl)-g_prime-P.delta_psl - P.gamma_sl*D2*Lsl;
    DODLlv = g_prime-P.delta_plv - P.gamma_lv*D2*Llv;
    f(:,2) = (-P.ksl)*DODLsl;
    f(:,1) = (D1*((h.^3/(3*P.eta)).*(D1*DODLlv)) + P.ksl*(P.delta_rho/P.rhol)*DODLsl);
return;

