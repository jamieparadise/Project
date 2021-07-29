function  plot_g(P)
    
        h=2.5/P.k1:0.1/P.k1:10/P.k1;
        g=EQ_g(h,P);
        plot(h*P.k1,g/P.C2,'LineWidth',1.5)    
        xlabel('$k_{1}h$','FontSize',18,'interpreter','latex')
        ylabel('$${C_2}^{-1}{g}$','FontSize',18,'interpreter','latex')
end

