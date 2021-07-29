function D2 = Centred_2nd_Deriv_5_Pt_Sym_Periodic(Nx,dx)
% Set up the second order differentiation matrix using central
% differences, with periodic boundary conditions.
% Create finite difference differentiation matrix
D2 = full(gallery('toeppen',Nx,-1,16,-30,16,-1));
% and impose BCs (they are periodic)
D2(Nx,1)=16;D2(Nx,2)=-1; 
D2(Nx-1,1)=-1;
D2(1,Nx)=16; D2(2,Nx)=-1;
D2(1,Nx-1)=-1;
D2 = D2/(12*dx^2); D2 = sparse(D2);
end

