classdef Simulation
    %SIMULATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Nx
        xL
        dx
        x
        dt
        L % layers at current time
        t=0; %current time
        L_history
        t_history
        D1
        D2
                
        pv %pressure in Pa
        T %temp in K
        Tc %temp in C
        %Parameter values  from Sup table 1
        C1 = 3.143e-3;
        C2 =4.116e-2;
        k1 = 2.043e9;
        k2 
        q0 = 7.148e9;
        alpha = 1.114;  %5.114;
        B = 7.785e-31;
        f= 1.106;
        a = 3.03e7;
        b= 5.0e8;
        
        %Parameter values  from Sup table 2
        Tt = 273.16;
        pt = 611.65;%
        delta_Hsv= 51059;%
        delta_Hlv = 45051;%
        delta_Hsl = 6008;
        
        %Parameter values  from Sup table 3
        R=8.31446261815324;
        dB = 0.37e-9;
        plv=453; %for the simulations in the paper
        

        rhol
        rhos 
        rhov 
        gamma_sl
        gamma_lv
        eta 
        rholv 
        u 
        klv %(1e-8)*S.rholv*S.T^(-1/2);
        ksl 
        delta_psl
        delta_plv
        delta_pk
        tau
        % More Parameters
        delta_rho
        qz
    end
    
    methods
        function S = Simulation(pv,T,Nx,dt,L0)
            S=Load_Physical_Parameters(S,pv,T);
            
            S.Nx = Nx; S.xL = 2500/S.k1;
            % dx is the the stepsize and x is the spatial domain
            S.dx = S.xL/(S.Nx-1); S.x = 0:S.dx:S.xL;
            S.dt=dt;
            
            % Create finite difference differentiation matrices (using central difference method), with periodic boundary conditions.
            %S.D1 = Central_Dif_1_Periodic(S.Nx,S.dx);
            %D1= Pt_5_Sym_Cent_Dif_1st_Periodic(Nx,dx);
            %S.D2 = Central_Dif_2_Periodic(S.Nx,S.dx);
            %D2= Pt_5_Sym_Cent_Dif_2nd_Periodic(Nx,dx);
            
            %Load the initial surface height of ice and water layers at the start time
            %Lsl and Llv are just the temporary surface heights at each time step
            S.L = L0;
            S.L_history(1,:)=S.L; 
            S.t_history(1)=0;
        end
        

        
        function S = Load_Physical_Parameters(S,pv,T)
            S.pv=pv; %pressure in Pa
            S.T=T; %temp in K
            S.Tc=S.T - 273.15; %temp in C
            
            %table 1
            S.k2 =2*S.k1;
            %table 3
            S.rhol = 55502 + 3.4549*S.Tc - 0.44461*S.Tc^2 + 0.0028885*S.Tc^3 -0.00031898*S.Tc^4;%
            S.rhos = 50885 - 9.71*S.Tc - 0.03*S.Tc^2;%
            S.rhov = S.pv/(S.R*S.T); % Rc or R?
            S.gamma_sl=(28 + 0.25*S.Tc)*10^-3; % 
            S.gamma_lv=(75.7 - 0.1775*S.Tc)*10^-3; % 
            S.eta = (1.39e-4)*(S.T/225-1)^-1.64; %
            S.rholv = S.plv/(S.R*S.T); % Rc or R?
            S.klv =(3.4e-10)*S.rhov*S.T^(-1/2); %(1e-8)*S.rholv*S.T^(-1/2)*10^-3;
            S.ksl=6.4*S.klv;
            S.u = ((S.dB/(2*pi))*S.rhos*S.delta_Hsl*2)/S.Tt;%1.3e-4
            %Parameter values  from Sup Note 6

            S.delta_psl=((S.rhos-S.rhol)*S.R*S.T*log(S.pv/S.pt)+(S.rhol*S.delta_Hlv-S.rhos*S.delta_Hsv)*(S.T-S.Tt)/S.Tt);%
            S.delta_plv=(S.rhol*S.R*S.T*log(S.pv/S.pt)-S.rhol*S.delta_Hlv*(S.T-S.Tt)/S.Tt);%
            %S.delta_pk=(S.rhos*S.ksl*S.delta_psl - S.rhol*S.klv*S.delta_plv)/(S.rhos*S.ksl+S.rhol*S.klv);
            S.tau = 3.8136e-11;% (3*P.eta)/(P.k1*P.gamma_lv)=3.8136e-11;%  using this formula. supposed to be 0.11e-9 in the paper
            % More parameters
            S.delta_rho = S.rhos-S.rhol;
            S.qz=2*pi/S.dB;
        end
        
        function S = Run_Simulation(run_time,time_method,odefun)
           S=Simulation;
           S.L=time_method(odefun,S.dt,round(run_time*S.tau/S.dt),S.L,S,S.D1,S.D2);
           S.t = S.t + run_time;
           S.L_history(end+1,:) = S.L;
           S.t_history(end+1) = S.t;
        end
        end
end

