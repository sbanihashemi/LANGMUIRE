%% this file generates input files for NHWAVE/SWAN 
clear;clc
% set cell numbers (one less than grid points)
Mglob = 6;
Nglob = 64;
Kglob = 60;

% cell size (in meter)
dx = 15.7067;
dy = 0.9817;

%% NHWAVE input file
% Both analytical initial condition and grid can be generated within NHWAVE
% but we can also use intial file and bathymetry file as input 

% initial file: eta0.txt
% NHWAVE reads file by line and put into x-direction first
Eta = zeros(Nglob,Mglob);
%save eta0.txt Eta -ASCII
%return
% initial file: uvw0.txt
% NHWAVE reads x-direction first, then y-direction, then z-direction
% for i = 1:Mglob;
%     for j=1:Nglob;
%         %sin(theta)
%         %Us(j,i,:) = u + sin((j-0.5)/Nglob*pi)*0.025*ones(Kglob,1);
%         %random noise
%         Us(j,i,:) = u + rand(Kglob,1)*0;
%     end
% end
% %return
% for k=1:Kglob;
%     U(((k-1)*Nglob+1):k*Nglob,:)=squeeze(Us(:,:,k));
% end
U = zeros(Nglob*Kglob,Mglob);
V = zeros(Nglob*Kglob,Mglob);
W = zeros(Nglob*Kglob,Mglob);
%save uvw0.txt U V W -ASCII 
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
%save swan_bathy.bot H_swan -ASCII
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
%save swan_grid_coord.grd swan_grid -ASCII



