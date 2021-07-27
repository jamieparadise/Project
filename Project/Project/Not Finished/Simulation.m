classdef Simulation
    %SIMULATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Nx
        xL
        dx
        x
        dt
        L
        
        D1
        D2
        time_method
        
        pv %pressure in Pa
        T %temp in K
        Tc %temp in C
        %Parameter values  from Sup table 1
        C1 = 3.143e-3
        C2 =4.116e-2
        k1 = 2.043e9
        k2 
        q0 = 7.148e9
        alpha = 1.514  %5.114;
        B = 7.785e-31
        f= 1.106
        a = 3.03e7
        b= 5.0e8
        
        %Parameter values  from Sup table 2
        Tt = 273.16;%
        %Ttc = 0.1;
        pt = 611.65;%
        % rhol = 55498;
        % rhos = 50888;
        % rhov  = 0.2694;
        delta_Hsv= 51059;%
        delta_Hlv = 45051;%
        %S.delta_Hsl = 6008;
        
        %Parameter values  from Sup table 3
        R=8.31446261815324
        dB = 0.37e-9%
        %Rc= 461.52;
        plv=517.5
        

        rhol
        rhos 
        rhov 
        gamma_sl
        gamma_lv
        eta 
        rholv 
        u = 1.3e-4%
        klv %(1e-8)*S.rholv*S.T^(-1/2);
        ksl 
        % More Parameters
        qz
    end
    
    methods
        function S = Simulation(pv,T,Nx,dt)
            S.Load_Physical_Parameters(pv,T);
            
            S.Nx = Nx; S.xL = 2500/S.k1;
            % dx is the the stepsize and x is the spatial domain
            S.dx = S.xL/(S.Nx-1); S.x = 0:S.dx:S.xL;
            S.dt=dt;
            
            % Create finite difference differentiation matrices (using central difference method), with periodic boundary conditions.
            S.D1 = Central_Dif_1_Periodic(S.Nx,S.dx);
            %D1= Pt_5_Sym_Cent_Dif_1st_Periodic(Nx,dx);
            S.D2 = Central_Dif_2_Periodic(S.Nx,S.dx);
            %D2= Pt_5_Sym_Cent_Dif_2nd_Periodic(Nx,dx);
            
            %Load the initial surface height of ice and water layers at the start time
            %Lsl and Llv are just the temporary surface heights at each time step
            S.L = INIT_Water_Drop_Ice_Flat(S.x,S);
            %L = INIT_Water_Flat_Ice_Square(x,P);
            
        end
        

        
        function S = Load_Physical_Parameters(pv,T)
            S.pv=pv; %pressure in Pa
            S.T=T; %temp in K
            S.Tc=S.T - 273.15; %temp in C
            
            %tabl 1
            S.k2 =2*S.k1;
            %table 3
            S.rhol = 55502 + 3.4549*S.Tc - 0.44461*S.Tc^2 + 0.0028885*S.Tc^3 -0.00031898*S.Tc^4;%
            S.rhos = 50885 - 9.71*S.Tc - 0.03*S.Tc^2;%
            S.rhov = S.pv/(S.R*S.T); % Rc or R?
            S.gamma_sl=(28 + 0.25*S.Tc)*10^-3; % 
            S.gamma_lv=(75.7 - 0.1775*S.Tc)*10^-3; % 
            S.eta = (1.39e-4)*(S.T/225-1)^-1.64; %
            S.rholv = S.plv/(S.R*S.T); % Rc or R?
            S.klv =(3.4e-10)*S.rhov*S.T^(-1/2)*10^-4; %(1e-8)*S.rholv*S.T^(-1/2);

            %Parameter values  from Sup Note 6

            S.delta_psl=((S.rhos-S.rhol)*S.R*S.T*log(S.pv/S.pt)+(S.rhol*S.delta_Hlv-S.rhos*S.delta_Hsv)*(S.T-S.Tt)/S.Tt);%
            S.delta_plv=(S.rhol*S.R*S.T*log(S.pv/S.pt)-S.rhol*S.delta_Hlv*(S.T-S.Tt)/S.Tt);%
            S.delta_pk=(S.rhos*S.ksl*S.delta_psl - S.rhol*S.klv*S.delta_plv)/(S.rhos*S.ksl+S.rhol*S.klv);
            S.tau = 0.11e-9;%(3*S.eta)/(S.k1*S.gamma_lv) = 3.8136e-11 using this formula gives 3.8136e-11
            % More parameters
            S.delta_rho = S.rhos-S.rhol;
            S.qz=2*pi/S.dB;
        end
    end
end

