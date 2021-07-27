function f = PDE_Full(t,L,P,D1,D2)
%formula from sup notes
    Nx=length(L)/2;
    Llv=L(1:Nx,1);
    Lsl=L(Nx+1:2*Nx,1);
    h=(Llv-Lsl); % h is also temporary
    gsr_prime =  -P.C2*P.k2*exp(-P.k2*h) ...
        +P.C1*exp(-P.k1*h).*(P.q0*sin(P.q0*h + P.alpha) + P.k1*cos(P.q0*h + P.alpha));
    gvdw_prime =  3*P.B*(h.^-4).*(1 - P.f*exp(-P.a*h)-(1-P.f)*exp(-P.b*h))...
        - P.B*(h.^-3).*( P.a*P.f*exp(-P.a*h)+P.b*(1-P.f)*exp(-P.b*h));
    g_prime = gsr_prime +gvdw_prime;
    DODLsl =  P.u*P.qz*sin(P.qz*Lsl)-P.delta_psl -g_prime-P.gamma_sl*D2*Lsl;
    DODLlv = g_prime-P.delta_plv - P.gamma_lv*D2*Llv;
    f(1:Nx,1) = (D1*((h.^3/(3*P.eta)).*(D1*DODLlv))- P.klv*DODLlv + P.ksl*(P.delta_rho/P.rhol)*DODLsl);
    f(Nx+1:2*Nx,1) = (-P.ksl)*DODLsl;
    return;
    
%     Llv=L(:,1);
%     Lsl=L(:,2);
%     h=(Llv-Lsl); % h is also temporary
%     gsr_prime =  -P.C2*P.k2*exp(-P.k2*h) ...
%         +P.C1*exp(-P.k1*h).*(P.q0*sin(P.q0*h + P.alpha) + P.k1*cos(P.q0*h + P.alpha));
%     gvdw_prime =  3*P.B*(h.^-4).*(1 - P.f*exp(-P.a*h)-(1-P.f)*exp(-P.b*h))...
%         - P.B*(h.^-3).*( P.a*P.f*exp(-P.a*h)+P.b*(1-P.f)*exp(-P.b*h));
%     g_prime = gsr_prime +gvdw_prime;
%     DODLsl = P.u*P.qz*sin(P.qz*Lsl)-P.delta_psl -g_prime- P.gamma_sl*D2*Lsl;
%     DODLlv = g_prime-P.delta_plv - P.gamma_lv*D2*Llv;
%     f(:,1) = (D1*((h.^3/(3*P.eta)).*(D1*DODLlv))- P.klv*DODLlv + P.ksl*(P.delta_rho/P.rhol)*DODLsl);
%     f(:,2) = (-P.ksl)*DODLsl;
% return;
