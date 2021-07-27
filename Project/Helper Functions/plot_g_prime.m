function  plot_g_prime(P)
    
        h=1.5/P.k1:0.1/P.k1:10/P.k1;
        g_prime=EQ_g_prime(h,P);
        plot(h*P.k1,g_prime/P.C2,'LineWidth',1.5)    
        xlabel('$k_{1}h$','FontSize',18,'interpreter','latex')
        ylabel('${C_2}^{-1}{g_prime}$','FontSize',18,'interpreter','latex')
end

