function L = multiplot(time_method,odefun,plot_times,dt,x,L,P,D1,D2)
%multiplot plots a figure made up of subplots in the same format as the plot_times matrix
total_run_time = 0;
[height,width]=size(plot_times);
for i = 1:width
    for j = 1:height
        if plot_times(j,i)>=0
            %run up to the next time value and plot the layers
            run_time = plot_times(j,i)-total_run_time;
            total_run_time= plot_times(j,i);
            L=time_method(odefun,dt,round(run_time*P.tau/dt),L,P,D1,D2);
            subplot(height,width,i+(j-1)*width);
            plot_layers(x,L,P);
            title("$t/tau=$" + plot_times(j,i),'interpreter','latex');
        elseif plot_times(j,i) == -1
            %skip to next plot
            continue
        elseif plot_times(j,i) == -2
            %plot some of the components of the DODLsl at the current time
            Nx= length(L)/2;
            subplot(height,width,i+(j-1)*width)
            plot(x*P.k1,ones(length(x))*(P.u*P.qz-P.delta_psl),'Color','r','LineStyle','--')
            hold on
            h = L(1:Nx) - L(Nx+1:end);
            plot(x*P.k1,(P.u*P.qz-P.delta_psl - EQ_g_prime(h,P)),'Color','r','LineStyle','-')
            hold on
            plot(x*P.k1,P.u*P.qz-P.delta_psl - EQ_g_prime(h,P) - P.gamma_sl*D2*L(Nx+1:end),'Color','black')
            
        end
    end
end
end

