function D1 = Centred_1st_Deriv_3_Pt_Sym_Periodic(Nx,dx)
% Set up the first order differentiation matrix using central
% differences, with periodic boundary conditions.
% Create finite difference differentiation matrix
D1 = full(gallery('tridiag',Nx,-1,0,1));
% and impose BCs (they are periodic)
D1(Nx,1)=1; D1(1,Nx)=-1;
D1 = D1/(2*dx); D1 = sparse(D1);
end

