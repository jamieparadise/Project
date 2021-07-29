function D1 = Centred_1st_Deriv_5_Pt_Sym_Periodic(Nx,dx)
% Set up the first order differentiation matrix using central
% differences, with periodic boundary conditions.
% Create finite difference differentiation matrix
D1 = full(gallery('toeppen',Nx,1,-8,0,8,-1));
% and impose BCs (they are periodic)
D1(Nx,1)=8;D1(Nx,2)=-1; 
D1(Nx-1,1)=-1;
D1(1,Nx)=-8; D1(2,Nx)=1;
D1(1,Nx-1)=1;
D1 = D1/(12*dx); D1 = sparse(D1);
end

