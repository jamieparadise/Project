function g = EQ_g(h,P)

    gsr =P.C2*exp(-P.k2*h) - P.C1*exp(-P.k1*h).*cos(P.q0*h + P.alpha);
    gvdw =  -P.B*(h.^-3).*(1 - P.f*exp(-P.a*h)-(1-P.f)*exp(-P.b*h));
    g = gsr +gvdw;
return;
