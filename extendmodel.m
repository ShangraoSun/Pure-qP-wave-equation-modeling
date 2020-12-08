% Extend Model
function [ext_model]=extendmodel(init_model,Nz,Nx,boundary_z,boundary_x)
Nx0 = Nx + 2*boundary_x;
Nz0 = Nz + 2*boundary_z;
ext_model = zeros(Nz0,Nx0);
% Region 1
for iz = boundary_z:Nz + boundary_z-1
    for ix = 1:boundary_x
        ext_model(iz,ix) = init_model(iz - boundary_z+1,1);
    end
end
% Region 2
for iz = boundary_z:Nz + boundary_z-1
    for ix = Nx + boundary_x:Nx + 2*boundary_x
        ext_model(iz,ix) = init_model(iz - boundary_z+1,Nx - 1);
    end
end
% Region 3
for iz = Nz + boundary_z:Nz + 2*boundary_z
    for ix = boundary_x:Nx + boundary_x-1
        ext_model(iz,ix) = init_model(Nz - 1,ix - boundary_x+1);
    end
end
% Region 4
for iz = 1:boundary_z
    for ix = boundary_x:Nx + boundary_x-1
        ext_model(iz,ix) = init_model(1,ix - boundary_x+1);
    end
end
% Region 5
for iz = 1:boundary_z
    for ix = 1:boundary_x
        ext_model(iz,ix) = init_model(1,1);
    end
end
% Region 6
for iz = 1:boundary_z
    for ix = Nx + boundary_x:Nx + 2*boundary_x
        ext_model(iz,ix) = init_model(1,Nx - 1);
    end
end
% Region 7
for iz = Nz + boundary_z:Nz + 2*boundary_z
    for ix = 1:boundary_x
        ext_model(iz,ix) = init_model(Nz - 1,1);
    end
end
% Region 8
for iz = Nz + boundary_z:Nz + 2*boundary_z
    for ix = Nx + boundary_x:Nx + 2*boundary_x
        ext_model(iz,ix) = init_model(Nz - 1,Nx - 1);
    end
end
% Region 9
for iz = boundary_z:Nz + boundary_z-1
    for ix = boundary_x:Nx + boundary_x-1
        ext_model(iz,ix) = init_model(iz - boundary_z+1, ix - boundary_x+1);
    end
end