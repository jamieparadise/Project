function  plot_layers(x,L,P)
            Nx= length(x);
            min_height_ice = min(L(Nx+1:end));
            first_lattice_line = floor(min_height_ice/P.dB);

            for k=1:30
               %This loop plots the red dashed lattice lines
               height_lattice_lines = (k-6+first_lattice_line)*P.k1*P.dB;
               plot(x*P.k1,ones(length(x),1)*height_lattice_lines,'Color','r','LineStyle','--','LineWidth',0.5)
               hold on
            end
            plot(x(:)*P.k1,L(Nx+1:end,end)*P.k1,'LineWidth',1.5,'Color','r');
            hold on
            plot(x(:)*P.k1,L(1:Nx,end)*P.k1,'LineWidth',1.5,'Color','b');
            xlabel('$k_{1}x$','FontSize',18,'interpreter','latex')
            ylabel('$k_{1}L$','FontSize',18,'interpreter','latex')
end