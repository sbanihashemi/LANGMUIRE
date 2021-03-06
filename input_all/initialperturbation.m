%% this file generates initial perturbation to NHWAVE
% domain size
% increase Mglob from 6 to 8 for SWAN MPI run 
clear;clc

Mglob = 4;
Nglob = 64;
Kglob = 60;

len_x = 94.20;
len_y = 62.83;
len_z = 15;

dx = len_x/Mglob;
dy = len_y/Nglob;


% eta0.txt 
% NHWAVE reads x-direction first
% No perturbation in eta0

eta0 = zeros(Nglob,Mglob);
save eta0.txt eta0 -ASCII

% uvw0.txt
load uprofstdy.mat
for i = 1:Mglob;
    for j=1:Nglob;
        %sin(theta)
        Us(j,i,:) = u + sin((j-0.5)/Nglob*pi)*0.025*ones(Kglob,1);
        %random noise
        %Us(j,i,:) = u + rand(Kglob,1)*0.025;
    end
end
%return
for k=1:Kglob;
    U(((k-1)*Nglob+1):k*Nglob,:)=squeeze(Us(:,:,k));
end
%V = rand(Nglob*Kglob,Mglob)*0.025;
% perturbation in V will shift the langmuir cells spanwise
% it is better not to use it 

V = zeros(Nglob*Kglob,Mglob);
W = zeros(Nglob*Kglob,Mglob);

save uvw0.txt U V W -ASCII 

%return
% bathymetry file: depth.txt
 for i = 1:Mglob;
         for j = 1:Nglob;
            H(j,i) = 15.0;
%           H(j,i) = 11.99 - 0.0125*20*(i-1);
         end
 end
%save depth.txt H -ASCII
%return
%% SWAN input file 
% SWAN input file is grid and bathymetry 
% bathymetry file: swan_bathy.bot
H_swan = H;
H_swan(H_swan<0)=0.1;
save swan_bathy.bot H_swan -ASCII
%return

% Generate SWAN computational grid and bathymetry
gridx = [0:Mglob-1]*dx;
gridy = [0:Nglob-1]*dy;
for j = 1:Nglob;
swan_gridx(((j-1)*Mglob+1):Mglob*j)=gridx;
end
MNglob = Mglob*Nglob;
for j = 1:Nglob;
    swan_gridy(((j-1)*Mglob+1):Mglob*j) = gridy(j)*ones(1,Mglob);
end
swan_grid(1:MNglob) = swan_gridx;
swan_grid(MNglob+1:2*MNglob) = swan_gridy;
swan_grid = swan_grid';
save swan_grid_coord.grd swan_grid -ASCII

