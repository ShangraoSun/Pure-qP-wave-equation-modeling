function [Dx2, Dz2, Dx1z1, Dx3z1, Dx1z3, Dx2z2, Dx4, Dz4] = FD_stencils(ix,iz,u,N,dx,dz)

% Part of code, just for finite difference (FD) stencils
% High-order spatial derivatives  and FD stencils

% dx, dz - Spatial step size
% Nx, Nz - Grid points
% N - Sptial order (twice of N - 2*N order)
% c1n - FD coefficients of spatial 1st-order derivative 
% c2n - FD coefficients of spatial 2nd-order derivative 

if N == 1   
    c1n = [0.500000000, 0.0, 0.0, 0.0, 0.0];
    c2n = [2.000000000, 0.0, 0.0, 0.0, 0.0];
elseif N == 2
    c1n = [0.666666667, -0.083333333, 0.0, 0.0, 0.0];
    c2n = [2.666666667, -0.166666667, 0.0, 0.0, 0.0];
elseif N == 3
    c1n = [0.750000000, -0.150000000, 0.016666667, 0.0, 0.0];
    c2n = [3.000000000, -0.300000000, 0.022222222, 0.0, 0.0];
elseif N == 4
    c1n = [0.800000000, -0.200000000, 0.038095238, -0.003571429, 0.0];
    c2n = [3.200000000, -0.400000000, 0.050793651, -0.003571429, 0.0];
elseif N == 5
    c1n = [0.833333333, -0.238095238, 0.059523810, -0.009920635, 0.000793651];
    c2n = [3.333333333, -0.476190476, 0.079365079, -0.009920635, 0.000634921];
end
% Coefficients of finite difference stencils of 1st-order spatial derivative
% 	c1n = [0.500000000, 0.0, 0.0, 0.0, 0.0];
% 	c1n = [0.666666667, -0.083333333, 0.0, 0.0, 0.0];
% 	c1n = [0.750000000, -0.150000000, 0.016666667, 0.0, 0.0];
% 	c1n = [0.800000000, -0.200000000, 0.038095238, -0.003571429, 0.0];
% 	c1n = [0.833333333, -0.238095238, 0.059523810, -0.009920635, 0.000793651];
% Coefficients of finite difference stencils of 2nd-order spatial derivative
% 	c2n = [2.000000000, 0.0, 0.0, 0.0, 0.0];
% 	c2n = [2.666666667, -0.166666667, 0.0, 0.0, 0.0];
% 	c2n = [3.000000000, -0.300000000, 0.022222222, 0.0, 0.0];
% 	c2n = [3.200000000, -0.400000000, 0.050793651, -0.003571429, 0.0];
% 	c2n = [3.333333333, -0.476190476, 0.079365079, -0.009920635, 0.000634921];

% The sum of c1n[ ] and c2n[ ]
c1n0 = 0.0; c2n0 = 0.0;
for i = 1:N
	c1n0 = c1n0 + c1n(i);
	c2n0 = c2n0 + c2n(i);
end

a1 = 0.0;
a2 = 0.0;
a3 = 0.0;
% Finite difference stencil for x's 2nd-order spatial derivative of wave field 'u'
% for ix = 1 + N:Nx - N
%     for iz = 1 + N:Nz - N
        for nx = 1:N
            a1 = a1 + c2n(nx)*(u(ix+nx,iz) + u(ix-nx,iz));
        end
        Dx2 = (c2n0*u(ix,iz) + a1)/(2*dx*dx);
%     end
% end
% Finite difference stencil for z's 2nd-order spatial derivative of wave field 'u'
% for ix = 1 + N:Nx - N 
%     for iz = 1 + N:Nz - N 
        for nz = 1:N
            a1 = a1 + c2n(nx)*(u(ix,iz+nz) + u(ix,iz-nz));
        end
        Dz2 = (c2n0*u(ix,iz) + a1)/(2*dz*dz);
%     end
% end
% Finite difference stencil for x's 1st-order and z's 1st-order spatial derivative of wave field 'u'
% for ix = 1 + N:Nx - N 
%     for iz = 1 + N:Nz - N 
        for nx = 1:N
            for nz = 1:N
                a1 = a1 + c1n(nx)*c1n(nz)*(u(ix+nx,iz+nz) - u(ix-nx,iz+nz) - u(ix+nx,iz-nz) - u(ix-nx,iz-nz));
            end
        end 
        Dx1z1 = a1/(dx*dz);
%     end
% end
% Finite difference stencil for x's 3th-order z's and 1st-order spatial derivative of wave field 'u'
Dx3z1 = ( c2n0*( c1n(1)*( c1n(1)*( u(iz+1,ix+1)- u(iz+1,ix-1)-u(iz-1,ix+1)+u(iz-1,ix-1) ) + c1n(2)*( u(iz+2,ix+1)- u(iz+2,ix-1)-u(iz-2,ix+1)+u(iz-2,ix-1) ) ) + ...
                            c1n(2)*( c1n(1)*( u(iz+1,ix+2)- u(iz+1,ix-2)-u(iz-1,ix+2)+u(iz-1,ix-2) ) + c1n(2)*( u(iz+2,ix+2)- u(iz+2,ix-2)-u(iz-2,ix+2)+u(iz-2,ix-2) ) ) ) + ...
                            c2n(1)*( c1n(1)*( c1n(1)*( u(iz+1,ix+2)- u(iz+1,ix)-u(iz-1,ix+2)+u(iz-1,ix) + u(iz+1,ix)-u(iz+1,ix-2)-u(iz-1,ix)+u(iz-1,ix-2) ) + ...
                                                          c1n(2)*( u(iz+2,ix+2)- u(iz+2,ix)-u(iz-2,ix+2)+u(iz-2,ix) + u(iz+2,ix)-u(iz+2,ix-2)-u(iz-2,ix)+u(iz-2,ix-2) ) ) + ...
                                           c1n(2)*( c1n(1)*( u(iz+1,ix+3)- u(iz+1,ix-1)-u(iz-1,ix+3)+u(iz-1,ix-1) + u(iz+1,ix+1)-u(iz+1,ix-3)-u(iz-1,ix+1)+u(iz-1,ix-3) ) + ...
                                                          c1n(2)*( u(iz+2,ix+3)- u(iz+2,ix-1)-u(iz-2,ix+3)+u(iz-2,ix-1) + u(iz+2,ix+1)-u(iz+2,ix-3)-u(iz-2,ix+1)+u(iz-2,ix-3) ) ) ) + ...
                            c2n(2)*( c1n(1)*( c1n(1)*( u(iz+1,ix+3)- u(iz+1,ix+1)-u(iz-1,ix+3)+u(iz-1,ix+1) + u(iz+1,ix-1)-u(iz+1,ix-3)-u(iz-1,ix-1)+u(iz-1,ix-3) ) + ...
                                                          c1n(2)*( u(iz+2,ix+3)- u(iz+2,ix+1)-u(iz-2,ix+3)+u(iz-2,ix+1) + u(iz+2,ix-1)-u(iz+2,ix-3)-u(iz-2,ix-1)+u(iz-2,ix-3) ) ) + ...
                                           c1n(2)*( c1n(1)*( u(iz+1,ix+4)- u(iz+1,ix)-u(iz-1,ix+4)+u(iz-1,ix) + u(iz+1,ix)-u(iz+1,ix-4)-u(iz-1,ix)+u(iz-1,ix-4) ) + ...
                                                          c1n(2)*( u(iz+2,ix+4)- u(iz+2,ix)-u(iz-2,ix+4)+u(iz-2,ix) + u(iz+2,ix)-u(iz+2,ix-4)-u(iz-2,ix)+u(iz-2,ix-4) ) ) ) ...
                  )/(2*dx*dx*dx*dz);
% Finite difference stencil for x's 1st-order z's and 3th-order spatial derivative of wave field 'u'
Dx1z3 = ( c2n0*( c1n(1)*( c1n(1)*( u(iz+1,ix+1)- u(iz+1,ix-1)-u(iz-1,ix+1)+u(iz-1,ix-1) ) + c1n(2)*( u(iz+2,ix+1)- u(iz+2,ix-1)-u(iz-2,ix+1)+u(iz-2,ix-1) ) ) + ...
                            c1n(2)*( c1n(1)*( u(iz+1,ix+2)- u(iz+1,ix-2)-u(iz-1,ix+2)+u(iz-1,ix-2) ) + c1n(2)*( u(iz+2,ix+2)- u(iz+2,ix-2)-u(iz-2,ix+2)+u(iz-2,ix-2) ) ) ) + ...
                            c2n(1)*( c1n(1)*( c1n(1)*( u(iz+2,ix+1)- u(iz+2,ix-1)-u(iz,ix+1)+u(iz,ix-1) + u(iz,ix+1)-u(iz,ix-1)-u(iz-2,ix+1)+u(iz-2,ix-1) ) + ...
                                                          c1n(2)*( u(iz+3,ix+1)- u(iz+3,ix-1)-u(iz-1,ix+1)+u(iz-1,ix-1) + u(iz+1,ix+1)-u(iz+1,ix-1)-u(iz-3,ix+1)+u(iz-3,ix-1) ) ) + ...
                                           c1n(2)*( c1n(1)*( u(iz+2,ix+2)- u(iz+2,ix-2)-u(iz,ix+2)+u(iz,ix-2) + u(iz,ix+2)-u(iz,ix-2)-u(iz-2,ix+2)+u(iz-2,ix-2) ) + ...
                                                          c1n(2)*( u(iz+3,ix+2)- u(iz+3,ix-2)-u(iz-1,ix+2)+u(iz-1,ix-2) + u(iz+1,ix+2)-u(iz+1,ix-2)-u(iz-3,ix+2)+u(iz-3,ix-2) ) ) ) + ...          
                           c2n(2)*( c1n(1)*( c1n(1)*( u(iz+3,ix+1)- u(iz+3,ix-1)-u(iz+1,ix+1)+u(iz+1,ix-1) + u(iz-1,ix+1)-u(iz-1,ix-1)-u(iz-3,ix+1)+u(iz-3,ix-1) ) + ...
                                                          c1n(2)*( u(iz+4,ix+1)- u(iz+4,ix-1)-u(iz,ix+1)+u(iz,ix-1) + u(iz,ix+1)-u(iz,ix-1)-u(iz-4,ix+1)+u(iz-4,ix-1) ) ) + ...
                                           c1n(2)*( c1n(1)*( u(iz+3,ix+2)- u(iz+3,ix-2)-u(iz+1,ix+2)+u(iz+1,ix-2) + u(iz-1,ix+2)-u(iz-1,ix-2)-u(iz-3,ix+2)+u(iz-3,ix-2) ) + ...
                                                          c1n(2)*( u(iz+4,ix+2)- u(iz+4,ix-2)-u(iz,ix+2)+u(iz,ix-2) + u(iz,ix+2)-u(iz,ix-2)-u(iz-4,ix+2)+u(iz-4,ix-2) ) ) ) ...
                  )/(2*dx*dz*dz*dz);
% Finite difference stencil for x's 2nd-order and z's 2nd-order spatial derivative of wave field 'u'
Dx2z2 = ( c2n0*( c2n0*( u(iz,ix) + c2n(1)*( u(iz+1,ix)-u(iz-1,ix) ) + c2n(2)*( u(iz+2,ix)-u(iz-2,ix) ) ) ) + ...
                c2n(1)*( c2n0*( u(iz,ix+1)+u(iz,ix-1) ) + c2n(1)*( u(iz+1,ix+1)+u(iz+1,ix-1)+u(iz-1,ix+1)+u(iz-1,ix-1) ) + ...
                                                                                    c2n(2)*( u(iz+2,ix+1)+u(iz+2,ix-1)+u(iz-2,ix+1)+u(iz-2,ix-1) ) ) + ...
                c2n(2)*( c2n0*( u(iz,ix+2)+u(iz,ix+2) ) + c2n(1)*( u(iz+1,ix+2)+u(iz+1,ix-2)+u(iz-1,ix+2)+u(iz-1,ix-2) ) + ...
                                                                                    c2n(2)*( u(iz+2,ix+2)+u(iz+2,ix-2)+u(iz-2,ix+2)+u(iz-2,ix-2) ) ) )/(4*dx*dx*dz*dz);
% Finite difference stencil for x's 4th-order  spatial derivative of wave field 'u'
Dx4 = ( c2n0*( c2n0*( u(iz,ix) + c2n(1)*( u(iz,ix+1)-u(iz,ix-1) ) + c2n(2)*( u(iz,ix+2)-u(iz,ix-2) ) ) ) + ...
                c2n(1)*( c2n0*( u(iz,ix+1)+u(iz,ix-1) ) + c2n(1)*( u(iz,ix+2)+u(iz,ix)+u(iz,ix)+u(iz,ix-2) ) + ...
                                                                                    c2n(2)*( u(iz,ix+3)+u(iz,ix-1)+u(iz,ix+1)+u(iz,ix-3) ) ) + ...
                c2n(2)*( c2n0*( u(iz,ix+2)+u(iz,ix+2) ) + c2n(1)*( u(iz,ix+3)+u(iz,ix+1)+u(iz,ix-1)+u(iz,ix-3) ) + ...
                                                                                    c2n(2)*( u(iz,ix+4)+u(iz,ix)+u(iz,ix)+u(iz,ix-4) ) ) )/(4*dx*dx*dx*dx);
% Finite difference stencil for z's 4th-order  spatial derivative of wave field 'u'
Dz4 = ( c2n0*( c2n0*( u(iz,ix) + c2n(1)*( u(iz+1,ix)-u(iz-1,ix) ) + c2n(2)*( u(iz+2,ix)-u(iz-2,ix) ) ) ) + ...
                c2n(1)*( c2n0*( u(iz+1,ix)+u(iz-1,ix) ) + c2n(1)*( u(iz+2,ix)+u(iz,ix)+u(iz,ix)+u(iz-2,ix) ) + ...
                                                                                    c2n(2)*( u(iz+3,ix)+u(iz-1,ix)+u(iz+1,ix)+u(iz-3,ix) ) ) + ...
                c2n(2)*( c2n0*( u(iz+2,ix)+u(iz+2,ix) ) + c2n(1)*( u(iz+3,ix)+u(iz+1,ix)+u(iz-1,ix)+u(iz-3,ix) ) + ...
                                                                                    c2n(2)*( u(iz+4,ix)+u(iz,ix)+u(iz,ix)+u(iz-4,ix) ) ) )/(4*dz*dz*dz*dz);
