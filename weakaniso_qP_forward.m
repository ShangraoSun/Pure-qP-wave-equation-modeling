%************************************************************************%
%***(TTI-Medium) weak anisotropic qP wave euqtion numerical simulation***%
%************************************************************************%

% Main Function
clear ; 
% close all;

dx = 10;     % spatial step length
dz = 10;
fm = 35.0;     % main-frequency 
dt = 0.001;     % time interval

boundary_z = 1;
boundary_x = 1;

Nt = 101;
Nx0 = 101;
Nz0 = 101;
Nx = Nx0 + 2*boundary_x;
Nz = Nz0 + 2*boundary_z;

time = Nt - 1;     % Time for snapshot

nsx = round((Nx-1)/2);  % Location of Source
nsz = round((Nz-1)/2);

%% Make Model
rho_init = 2.16*ones(Nz0,Nx0);
vp0_init = 3410*ones(Nz0,Nx0);
vs0_init = 2070*ones(Nz0,Nx0);
epsilon_init = 0.2*ones(Nz0,Nx0);
delta_init = 0.2*ones(Nz0,Nx0);
theta0_init = 0*pi/180*ones(Nz0,Nx0);

%% Extend Model
[rho]=extendmodel(rho_init,Nz0,Nx0,boundary_z,boundary_x);
[vp0]=extendmodel(vp0_init,Nz0,Nx0,boundary_z,boundary_x);
[vs0]=extendmodel(vs0_init,Nz0,Nx0,boundary_z,boundary_x);
[epsilon]=extendmodel(epsilon_init,Nz0,Nx0,boundary_z,boundary_x);
[delta]=extendmodel(delta_init,Nz0,Nx0,boundary_z,boundary_x);
[theta0]=extendmodel(theta0_init,Nz0,Nx0,boundary_z,boundary_x);

%% Make Sparse Matrix
tic

c1 = 1/(dx^2);    % coefficient
c2 = 1/(dz^2);
c0 = c1 + c2;

n = Nx^2;
N = 1;

C1 = ones(n,1)*c1; 
C2 = ones(n,1)*c2;
C0 = -2*(C1+C2);

for i=2:n-2
    if mod(i,Nx) == 0
        C1(i) = 0;
    end
end

C00 = spdiags(C0, 0, n, n);
C10 = spdiags(C1, -1, n, n);
C20 = spdiags(C2, -Nx, n, n);
C = C00 + C10 + C20 + C10' + C20';

disp('***********Complete Matrix C Assignment************');
toc

%% Source 
source = zeros(Nt,1);
for it = 1:Nt
    t = (it - 1.2/fm/dt)*dt;
    source(it) = (1-2*(pi*fm*t)^2)*exp(-(pi*fm*t)^2);
end
disp('***********Complete The "Source" Assignment************');

%% File name
Snap_FileName = './qP_Datas/iso_qP_snap_Nx%d_Nz%d_Nt%d.dat';
Snap_FileName0 = sprintf(Snap_FileName, Nz,Nx, Nt);

RecordX_FileName = './qP_Datas/iso_qP_recordX_Nx%d_Nz%d_Nt%d.dat';
RecordX_FileName0 = sprintf(RecordX_FileName, Nz,Nx, Nt);

RecordZ_FileName = './qP_Datas/iso_qP_recordZ_Nx%d_Nz%d_Nt%d.dat';
RecordZ_FileName0 = sprintf(RecordZ_FileName, Nz,Nx, Nt);

snap_filename = './qP_Datas/iso_qP_snap_Nx%d_Nz%d_Nt%d.mat';
snap_filename0 = sprintf(snap_filename, Nz,Nx, Nt);

recordX_filename = './qP_Datas/iso_qP_recordX_Nx%d_Nz%d_Nt%d.mat';
recordX_filename0 = sprintf(recordX_filename, Nz,Nx, Nt);

recordZ_filename = './qP_Datas/iso_qP_recordZ_Nx%d_Nz%d_Nt%d.mat';
recordZ_filename0 = sprintf(recordZ_filename, Nz,Nx, Nt);

%% Begin to Calculate
Pp = zeros(Nz,Nx);
Qp = zeros(Nz,Nx);
Qpm = zeros(Nz,Nx);
Qpp = zeros(Nz,Nx);
record_qPX = zeros(Nt,Nx);
record_qPZ = zeros(Nt,Nz);
%% Absorbing Boundary Condition (ABCs, Cerjan,1998)
apara = 0.01;
ABCs = ones(Nz,Nx);

% Z-direction
for iz=1:Nz
    for ix=1:Nx
        if(iz < boundary_z)
            ABCs(iz,ix) = ABCs(iz,ix)*exp( -( ( apara * (boundary_z - iz) )^3 ) );
        elseif(iz > (Nz - boundary_z + 1))
            ABCs(iz,ix) = ABCs(iz,ix)*exp( -( ( apara * (iz - Nz + boundary_z - 1) )^3 ) );
        end
    end
end
% X-direction
for iz=1:Nz
    for ix=1:Nx
        if(ix < boundary_x)
            ABCs(iz,ix) = ABCs(iz,ix)*exp( -( ( apara * (boundary_x - ix) )^3 ) );
        elseif(ix > Nx - boundary_x + 1)
            ABCs(iz,ix) = ABCs(iz,ix)*exp( -( ( apara * (ix - Nx + boundary_x - 1) )^3 ) );
        end
    end
end

% return;
%%
tic
disp('*************** Ready to run !!!****************');
for it = 1:Nt
    
    if mod(it,50) == 0 || it == Nt 
        fprintf("Time = %d,\tTravel time: %f\n", it, it*dt);
    elseif it == Nt-1
        fprintf("Complete...")
    end
        
    for iz = 1+2:Nz - 2
        for ix = 1+2:Nx - 2
            
            Dx4 = (Pp(iz,ix+2)-4*Pp(iz,ix+1)+6*Pp(iz,ix)-4*Pp(iz,ix-1)+Pp(iz,ix-2))/(dx^4);  % 4-th Spatial difference
            Dz4 = (Pp(iz+2,ix)-4*Pp(iz+1,ix)+6*Pp(iz,ix)-4*Pp(iz-1,ix)+Pp(iz-2,ix))/(dz^4);
            Dxz2 = (4*Pp(iz,ix)-2*(Pp(iz,ix-1)+Pp(iz-1,ix)+Pp(iz,ix+1)+Pp(iz+1,ix))+Pp(iz-1,ix-1)+Pp(iz-1,ix+1)+Pp(iz+1,ix+1)+Pp(iz+1,ix-1))/dx^2/dz^2;
            Dx1z3 = (Pp(iz+2,ix+1)-Pp(iz-2,ix+1)-Pp(iz+2,ix-1)+Pp(iz-2,ix-1)+2*(Pp(iz+1,ix-1)-Pp(iz-1,ix-1)-Pp(iz+1,ix+1)+Pp(iz-1,ix+1)))/(4*dx*dz^3);       
            Dx3z1 = (Pp(iz+1,ix+2)-Pp(iz+1,ix-2)-Pp(iz-1,ix+2)+Pp(iz-1,ix-2)+2*(Pp(iz-1,ix+1)-Pp(iz-1,ix-1)-Pp(iz+1,ix+1)+Pp(iz+1,ix-1)))/(4*dx^3*dz);
            % Version 2  (TTI qP wave) - Liang Kai et al., 2009, Wave equation decomposition in 3-D TTI medium
            coeffiqP1 = 1+2*epsilon(iz,ix)*cos(theta0(iz,ix))^4+0.5*delta(iz,ix)*sin(2*theta0(iz,ix)^2);  %Dx4
            coeffiqP2 = 1+2*epsilon(iz,ix)*sin(theta0(iz,ix))^4+0.5*delta(iz,ix)*sin(2*theta0(iz,ix)^2);  %Dz4
            coeffiqP3 = 2*(1+delta(iz,ix))+3*(epsilon(iz,ix)-delta(iz,ix))*sin(2*theta0(iz,ix)^2); %Dxz2
            coeffiqP4 = 8*epsilon(iz,ix)*cos(theta0(iz,ix))^3*sin(theta0(iz,ix))-delta(iz,ix)*sin(4*theta0(iz,ix));   %Dx3z1
            coeffiqP5 = 8*epsilon(iz,ix)*cos(theta0(iz,ix))*sin(theta0(iz,ix))^3+delta(iz,ix)*sin(4*theta0(iz,ix));   %Dx3z1      
            Qpp(iz,ix) = 2*Qp(iz,ix) - Qpm(iz,ix) + dt^2*vp0(iz,ix)^2*...
                        ( coeffiqP1*Dx4+coeffiqP2*Dz4+coeffiqP3*Dxz2+coeffiqP4*Dx3z1+coeffiqP5*Dx1z3 );
            
        end
    end
    
    for iz = 1:Nz
        for ix = 1:Nx
            Qpm(iz,ix) = Qp(iz,ix)*ABCs(iz,ix);		   % ABCs%
            Qp(iz,ix) = Qpp(iz,ix)*ABCs(iz,ix);
        end
    end
    
    Qp(nsz,nsx) = Qp(nsz,nsx) + source(it);    % Source Loading
    
    for ix = 1:Nx
        record_qPX(it,ix) = Qp(nsz,ix);
    end
    
    for iz = 1:Nz
        record_qPZ(it,iz) = Qp(iz,nsx);
    end
    
    % Calculate the wave field value
    Pp = reshape(C\(Qp(:)),[Nz,Nx]);   
    
     if it == time
        
        fids = fopen(Snap_FileName0,'wb');   % Snaoshot 
        fwrite(fids,Qp,'float');
        fclose(fids);

        fidr = fopen(RecordX_FileName0,'wb');   % X-direction Record 
        fwrite(fidr,record_qPX,'float');
        fclose(fidr);
        
        fidr = fopen(RecordZ_FileName0,'wb');   % Z-direction Record 
        fwrite(fidr,record_qPZ,'float');
        fclose(fidr);
        
     end
    
%     if rem(it,100) == 0
%         figure;
%         imagesc(Qpp(boundary_z:boundary_z+Nz0,boundary_x:boundary_x+Nx0));
%         axis square; colormap(gray);
%         set(gca,'fontsize',16);
%         set(gcf,'position',[200,200,500,500]); 
%     end
    
end  

toc 

save(snap_filename0,'Qp');

save(recordX_filename0,'record_qPX');

save(recordZ_filename0,'record_qPZ');


figure;
imagesc(Qpp(boundary_z:boundary_z+Nz0,boundary_x:boundary_x+Nx0));
axis square; colormap(gray);
set(gca,'fontsize',16);
set(gcf,'position',[200,200,500,500]); 

figure;
imagesc(record_qPX(:,boundary_x:boundary_x+Nx0));
axis square; colormap(gray);
set(gca,'fontsize',16);
set(gcf,'position',[200,200,500,500]); 

figure;
imagesc(record_qPZ(:,boundary_z:boundary_z+Nz0));
axis square; colormap(gray);
set(gca,'fontsize',16);
set(gcf,'position',[200,200,500,500]); 

disp('****************** Am done !!! ********************');

            