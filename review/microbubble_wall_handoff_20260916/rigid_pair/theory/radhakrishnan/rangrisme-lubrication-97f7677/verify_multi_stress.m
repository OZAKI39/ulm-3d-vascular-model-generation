%-------------------------------------------------------------------------%
% Script to verify stress contributions due to multiple particles
% Date May 1, 2017
% Modified by Ranga for the new lubrication class Simple
%-------------------------------------------------------------------------%

clear;
%close all;

display('Reading parameters');
%muf = nuf*rhof;
mu = 1; %FLD_VISC


Kn=10000;

folder = '/tmp/ranga/two_particles/Simple/';
filename2 = strcat(folder,'dump.lammpstrj');

gdot =0.0; % 0.05;

N = 51; %Number of steps

nPart = 2;


dt=1e-5;

%% Read in s.p

display('Reading particle output file');

f1 = 'timestep'; v1 = zeros(1,1);
f2 = 'natoms'; v2 = zeros(1,1);
f3 = 'xbox'; v3 = zeros(1,3);
f4 = 'ybox'; v4 = zeros(1,3);
f5 = 'zbox'; v5 = zeros(1,3);
f6 = 'values'; v6 = zeros(1,19);


pdata(1:N,1:nPart) = struct(f1, v1, f2, v2, f3, v3, f4, v4, f5, v5, f6, v6);


file = fopen(filename2);

for i=1:N
    pdata(i,1).timestep = fscanf(file, '%*s %*s\n%d\n', 1);
    pdata(i,1).natoms = fscanf(file, '%*s %*s %*s %*s\n%d\n', 1);
    pdata(i,1).xbox = fscanf(file, '%*s %*s %*s %*s %*s %*s %*s %*s %*s\n%f %f %f\n', [1, 3]);
    pdata(i,1).ybox = fscanf(file, '%f %f %f\n', [1, 3]);
    pdata(i,1).zbox = fscanf(file, '%f %f %f\n', [1, 3]);
    
    for j=1:pdata(i,1).natoms
        if j == 1
            %dump id2 all custom 100 s.p id type radius mass x y z vx vy vz fx fy fz tqx tqy tqz omegax omegay omegaz
            tmp = fscanf(file, '%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', [1 19]);
        else
            tmp = fscanf(file, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', [1 19]);
        end
        pdata(i,tmp(1)).values = tmp;
        
    end
    
end
fclose(file);


%% Some Variable declaration
display('Declaring variables');
fytot=zeros([N nPart]);
fylub=zeros([N nPart]);
fyfric=zeros([N nPart]);

xpos=zeros([N nPart]);
ypos=zeros([N nPart]);
zpos=zeros([N nPart]);

Ux=zeros([N nPart]);
Uy=zeros([N nPart]);
Uz=zeros([N nPart]);

Fx=zeros([N nPart]);
Fy=zeros([N nPart]);
Fz=zeros([N nPart]);
type=zeros([N nPart]);
time=zeros([N 1]);
f_dimensionless = zeros(N,1);
flubanalyticaln = zeros(N,1);

torx = zeros(N,1);
tory = zeros(N,1);
torz = zeros(N,1);

toranalx=zeros(N,1);
toranaly=zeros(N,1);
toranalz=zeros(N,1);

for i = 1:N %size(s,2)
    for j=1:pdata(i,1).natoms
        xpos(i,j)=pdata(i,j).values(5);
        ypos(i,j)=pdata(i,j).values(6);
        zpos(i,j)=pdata(i,j).values(7);
        
        type(i,j)=pdata(i,j).values(2);
        
        Ux(i,j)=pdata(i,j).values(8);
        Uy(i,j)=pdata(i,j).values(9);
        Uz(i,j)=pdata(i,j).values(10);
        
        Fx(i,j)=pdata(i,j).values(11);
        Fy(i,j)=pdata(i,j).values(12);
        Fz(i,j)=pdata(i,j).values(13);
    end
    time(i)=dt*pdata(i,1).timestep;
end



lx = pdata(1,1).xbox(2)- pdata(1,1).xbox(1);
ly = pdata(1,1).ybox(2)- pdata(1,1).ybox(1);
lz = pdata(1,1).zbox(2)- pdata(1,1).zbox(1);
box_vol=lx*ly*lz;

initial=dt*0;
delx = (time-initial)*gdot*ly; % assume shear rate = dvx/dy box displacement
%delx = disp; % -L/2 <delx < L/2
delv = gdot*ly;

Scontxy=zeros(N,1);
Slubxy=zeros(N,1);

%% Computation of force for comparison
display('Computing forces')

for i = 1:N %size(s,2)
    f_dimensionless(i) =sqrt(Fx(i,1).^2+Fy(i,1).^2+Fz(i,1).^2);
    %f_dimensionless(i) =Fy(i,1);
    torx(i)=pdata(i,1).values(14);
    tory(i)=pdata(i,1).values(15);
    torz(i)=pdata(i,1).values(16);
    for j=1:pdata(i,1).natoms
        r1=pdata(i,j).values(3);
        
        for k=j+1:pdata(i,1).natoms
            r2=pdata(i,k).values(3);
            %calculate contact forces
            
            [ccdx21,ccdy21,ccdz21,ccd21]=minimage(lx,ly,lz,delx(i),xpos(i,j),ypos(i,j),zpos(i,j),xpos(i,k),ypos(i,k),zpos(i,k));
            
            nx=ccdx21/ccd21;
            ny=ccdy21/ccd21;
            nz=ccdz21/ccd21;
            eps21 = ccd21 - r1 - r2;
            if eps21<0
                fx1=Kn*eps21*nx;
                fy1=Kn*eps21*ny;
                fz1=Kn*eps21*nz;
            else
                fx1=0;
                fy1=0;
                fz1=0;
            end
            
            %Sum stress contribution due to contact forces
            Scontxy(i)=Scontxy(i)+fy1*ccdx21;
            
            
            %Calculate the lubrication force
            %[fx2,fy2]=calc_lub();
            %r12= r1.*r2./(r1+r2);
            hsep=2*eps21/(r1+r2);
            hcutin=0.001*2;
            hcutout=0.05*2;
            beta=r2/r1;
            XA11_a=6*pi*r1*(2*beta^2/(1+beta)^3);
            XA11_b=6*pi*r1*(beta*(1+7*beta+beta^2)/(5*(1+beta)^3));
            YA11_b=6*pi*r1*(4*beta*(2+beta+2*beta^2)/(15*(1+beta)^3));
            YB11_b=-(4*pi*r1^2)*beta*(4+beta)/(5*(1+beta)^2);
            YC12_b=0.8*pi*r1^3*beta^2/(1+beta)*(1-4.0/beta);
            
            
            Omg1x= pdata(i,j).values(17);
            Omg1y= pdata(i,j).values(18);
            Omg1z= pdata(i,j).values(19);
            
            Omg2x= pdata(i,k).values(17);
            Omg2y= pdata(i,k).values(18);
            Omg2z= pdata(i,k).values(19);
            
            wsx=(Omg1x+Omg2x)/2.0;
            wsy=(Omg1y+Omg2y)/2.0;
            wsz=(Omg1z+Omg2z)/2.0;
            wsnr=(nx*wsx + ny*wsy+nz*wsz);
            
            wdx=(Omg1x-Omg2x)/2.0;
            wdy=(Omg1y-Omg2y)/2.0;
            wdz=(Omg1z-Omg2z)/2.0;
            wdnr=(nx*wdx + ny*wdy+nz*wdz);
            
            wsnz=wsx*ny-wsy*nx;
            wsnx=wsy*nz-wsz*ny;
            wsny=wsz*nx-wsx*nz;
           
            
            wdnz=wdx*ny-wdy*nx;
            wdnx=wdy*nz-wdz*ny;
            wdny=wdz*nx-wdx*nz;
            
            
            Vx=min_velx(ypos(i,j),ypos(i,k),Ux(i,j),Ux(i,k),ly,delv);
            Vy=Uy(i,k)-Uy(i,j);
            Vz=Uz(i,k)-Uz(i,j);
            
            
            Vnnr=(nx*Vx + ny*Vy+nz*Vz);
            Unx = Vnnr.*nx;
            Uny = Vnnr.*ny;
            Unz= Vnnr.*nz;
            
            Utx = Vx - Unx;
            Uty = Vy - Uny;
            Utz = Vz -Unz;
            
            
            if hsep<=hcutout
                if hsep<=hcutin
                    hsep=hcutin;
                end
                
                F_l_nx =  (XA11_a/hsep+XA11_b*log(1/hsep))*Unx;
                F_l_tx = YA11_b*Utx*log(1.0/hsep);
                F_l_wx1 =-(r1+r2)*YA11_b*wsnx*log(1.0/hsep);
                F_l_wx2 =+(1-beta*(1+4*beta)/(4+beta))*YB11_b*wdnx*log(1.0/hsep);
                
                F_l_ny =  (XA11_a/hsep+XA11_b*log(1/hsep))*Uny;
                F_l_ty = YA11_b*Uty*log(1.0/hsep);
                F_l_wy1 = -(r1+r2)*YA11_b*wsny*log(1.0/hsep);
                F_l_wy2 =+(1-beta*(1+4*beta)/(4+beta))*YB11_b*wdny*log(1.0/hsep);
                
                F_l_nz =  (XA11_a/hsep+XA11_b*log(1/hsep))*Unz;
                F_l_tz = YA11_b*Utz*log(1.0/hsep);
                F_l_wz1 = -(r1+r2)*YA11_b*wsnz*log(1.0/hsep);
                F_l_wz2 =+(1-beta*(1+4*beta)/(4+beta))*YB11_b*wdnz*log(1.0/hsep);
                
                T_l_nx=mu*(YB11_b*( Vy*nz -Vz*ny + (r1+r2)*(wsx-wsnr*nx) ) + YC12_b*(wdx-wdnr*nx) )*log(1/hsep);
                T_l_ny=mu*(YB11_b*( Vz*nx -Vx*nz + (r1+r2)*(wsy-wsnr*ny) ) + YC12_b*(wdy-wdnr*ny) )*log(1/hsep);
                T_l_nz=mu*(YB11_b*( Vx*ny -Vy*nx + (r1+r2)*(wsz-wsnr*nz) ) + YC12_b*(wdz-wdnr*nz) )*log(1/hsep);
                
            else
                F_l_nx=  0;
                F_l_tx = 0;
                F_l_wx1 =0;
                F_l_wx2 =0;
                F_l_ny = 0;
                F_l_ty = 0;
                F_l_wy1 =0;
                F_l_wy2 =0;
                T_l_nx=0;
                T_l_ny=0;
            end
            
            
            fx2=mu*(F_l_nx+F_l_tx+F_l_wx1+F_l_wx2);
            fy2=mu*(F_l_ny+F_l_ty+F_l_wy1+F_l_wy2);
            fz2=mu*(F_l_nz+F_l_tz+F_l_wz1+F_l_wz2);
            
            flubanalyticaln(i)=sqrt((fx1+fx2).^2+(fy1+fy2).^2+(fz1+fz2).^2);
            
            toranalx(i)=T_l_nx;
            toranaly(i)=T_l_ny;
            toranalz(i)=T_l_nz;
            
            Slubxy(i)=Slubxy(i)+fy2*ccdx21;
            
        end
        
    end
    
end


for i = 1:N %size(s,2)
    
    Scontxy(i)=Scontxy(i)/box_vol;
    Slubxy(i)=Slubxy(i)/box_vol;
end



%% Plotting of solution

display('Plot of analytical solution for stress')


    subplot(3,1,1);
   plot(time, f_dimensionless, 'o','LineWidth', 2);
    hold on;
    plot(time, flubanalyticaln, '-r','linewidth', 2);
    
    %xlabel('Time');
    ylabel('Force');
    set(gca,'XScale','log');
    %hold off;
    
    
    subplot(3,1,2);
    plot(time,sqrt(toranalx.^2+toranaly.^2+toranalz.^2), '-b','LineWidth', 2);
    hold on;
     plot(time,sqrt(torx.^2+tory.^2+torz.^2), 'xb','LineWidth', 2);
    %xlabel('Time');
    ylabel('Torque');
    set(gca,'XScale','log');
    %hold off;
    
    
    
    %plot(time,-SC_t, '-k','LineWidth', 2);
    
    
    %plot(time,-SC_con, '-r','LineWidth', 2);
    subplot(3,1,3);


sdata=dlmread(strcat(folder,'stressOutputKE'), ' ', 1, 0);
plot(time(2:end),Slubxy(2:end), '-b','LineWidth', 2);
hold all;
% plot(time(2:end),Scontxy(2:end), '-r','LineWidth', 2);
% plot(time(2:end),Slubxy(2:end)+Scontxy(2:end), '-k','LineWidth', 2);
% 
plot(sdata(:,1)*dt,-sdata(:,2), 'ok','LineWidth', 2);
% 
% plot(sdata(:,1)*dt,-sdata(:,3), 'or','LineWidth', 2);
% plot(sdata(:,1)*dt,-sdata(:,4), 'ob','LineWidth', 2);
% 
% 
set(gca,'XScale','log');
xlabel('Time');
ylabel('Stress');
%hold off;
