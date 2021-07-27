function  plot_wk(P)
    
        h=1.5/P.k1:0.1/P.k1:10/P.k1;
        wk=EQ_wk(h,P);
        plot(h*P.k1,wk/P.C2,'LineWidth',1.5)    
        xlabel('$k_{1}h$','FontSize',18,'interpreter','latex')
        ylabel('${C_2}^{-1}{\omega}_{k}$','FontSize',18,'interpreter','latex')
end

