function g_prime = EQ_g_prime(h,P)
    gsr_prime =  -P.C2*P.k2*exp(-P.k2*h) ...
        +P.C1*exp(-P.k1*h).*(P.q0*sin(P.q0*h + P.alpha) + P.k1*cos(P.q0*h + P.alpha));
    gvdw_prime =  3*P.B*(h.^-4).*(1 - P.f*exp(-P.a*h)-(1-P.f)*exp(-P.b*h))...
        - P.B*(h.^-3).*( P.a*P.f*exp(-P.a*h)+P.b*(1-P.f)*exp(-P.b*h));

    g_prime = gsr_prime +gvdw_prime;
return;
