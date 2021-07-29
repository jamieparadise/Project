function ICE_main
%add all subfolders to search path
addpath(genpath(pwd));

% P is the structure array containing physical parameters
% Load_Parameters takes pressure in Pa and temperature in K as arguments
P = Load_Parameters(518.5, 269.5);

% Nx is the number of spacial steps, xL the size of the domain
Nx = 250; xL = 2500/P.k1;
% dx is the the stepsize and x is the spatial domain
dx = xL/(Nx-1); x = 0:dx:xL;


% Load finite difference differentiation matrices (using central difference method), with periodic boundary conditions.
D1 = Centred_1st_Deriv_3_Pt_Sym_Periodic(Nx,dx);
%D1= Pt_5_Sym_Cent_Dif_1st_Periodic(Nx,dx);
D2 = Centred_2nd_Deriv_3_Pt_Sym_Periodic(Nx,dx);
%D2= Pt_5_Sym_Cent_Dif_2nd_Periodic(Nx,dx);

dt =1e-12;  P.dt =dt;                % timestep

%plot_g(P);
%Load the initial surface height of ice and water layers at the start time
%L(1,:)= Llv and L(2,:)=Lsl which are just the temporary surface heights at each time step
L0 = INIT_Water_Drop_Ice_Flat(x,P);
%L0 = INIT_Water_Flat_Ice_Square(x,P);

%S=Simulation(518.5,269.5,Nx,dt,L0)

%this runs up a time (in terms of t/tau, like in the paper) and plots the layers
%run_time = 3.3e4;
%L=Forward_RK4(@SDE_Full_w_spectral_d,dt,round(run_time*P.tau/dt),L0,P,D1,D2);
%L=Forward_RK4(@PDE_Full,dt,round(run_time*P.tau/dt),L0,P,D1,D2);

% plot_layers(x,L,P);

%This will run the simulation up to the last (positive) number in the plot_times matrix.
%It will also plot the layers at the (positive) times listed and plot components of DODLsl in figures to the right.  
%figure('position',[100 100 1600 500]);
plot_times=[1, 4.3e4,6.1e5;
            -2,-2,-2];
L= multiplot(@Forward_RK4,@SDE_Full_w_spectral_d,plot_times,dt,x,L0,P,D1,D2);

% plot_times=[1, 5e2,5e3;
%             -2,-2,-2;
%             -1,-1,-1;
%             -1,-1,-1];
% L= multiplot(@Forward_RK4,@SDE_Full_w_spectral_d,plot_times,dt,x,L0,P,D1,D2);
% plot_times=[-1,-1,-1;
%             -1,-1,-1;
%             1, 5e2,5e3;
%             -2,-2,-2];
% L= multiplot(@Forward_RK4,@PDE_Full_w_spectral_d,plot_times,dt,x,L0,P,D1,D2);
%mass=TEST_is_mass_conserved(0.01,@Euler_Method,@PDE_No_Evap_Or_Cond,10000,dt,x,L,P,D1,D2)
%save2pdf("Ice_PDE_(Nt="+Nt+")(dt="+dt+")",gcf,300);

end