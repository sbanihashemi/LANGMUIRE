clc;clear

size_x = 4;
size_y = 128;
size_z = 60;

len_x = 94.24*8/6;
len_y = 62.83;
len_z = 15;

dx= len_x/size_x;
dy= len_y/size_y;
dz= len_z/size_z;

x = [dx:dx:len_x];
y = [dy:dy:len_y];
z = -[len_z:-dz:dz];

time=load('time');
time = round(time);
t_str = int2str(time);

nFrames=length(time);
%return

% load data

for i=1:nFrames;
%for i=45
    
    fnum = sprintf('%.4d',i);
    uu=load(['u_' fnum]);
    vv=load(['v_' fnum]);
    ww=load(['w_' fnum]);
    eta=load(['eta_' fnum]);
    
    u_surf = uu((size_z-5-1)*size_y+1:(size_z-5)*size_y,:);
    v_surf = vv((size_z-5-1)*size_y+1:(size_z-5)*size_y,:);
    w_surf = ww((size_z-5-1)*size_y+1:(size_z-5)*size_y,:);
    
    u_mid = uu((size_z/2-1)*size_y+1:size_z/2*size_y,:);
    v_mid = vv((size_z/2-1)*size_y+1:size_z/2*size_y,:);
    w_mid = ww((size_z/2-1)*size_y+1:size_z/2*size_y,:);
    
    u_bot = uu((size_z/12-1)*size_y+1:size_z/12*size_y,:);
    v_bot = vv((size_z/12-1)*size_y+1:size_z/12*size_y,:);
    w_bot = ww((size_z/12-1)*size_y+1:size_z/12*size_y,:);
   
%end
    
    umean_surf = mean(u_surf(:,1));
    vmean_surf = mean(v_surf(:,1));
    
    udiff_surf = u_surf-umean_surf;
    vdiff_surf = v_surf-vmean_surf;
    
    udiff_min = min(min(udiff_surf));
    udiff_max = max(max(udiff_surf));
    
    vdiff_min = min(min(vdiff_surf));
    vdiff_max = max(max(vdiff_surf));
    
%     figure(1);
%     clf;
%     subplot(2,1,1)
%     contourf(x,y,udiff_surf);
%     %caxis([udiff_min udiff_max])
%     caxis([-0.05 0.05])
%     colorbar
%     xlabel('X (m)','FontSize',10);
%     ylabel('Y (m)','FontSize',10);
%     title(['U (Time = ',t_str(i,:),' s)'])
%     
%     subplot(2,1,2)
%     contourf(x,y,vdiff_surf);
%     caxis([-0.05 0.05])
%     %caxis([vdiff_min vdiff_max])
%     colorbar
%     xlabel('X (m)','FontSize',10);
%     ylabel('Y (m)','FontSize',10);
%     title(['V (Time = ',t_str(i,:),' s)'])
%     
%     pause(0.5)
% x-direction uniform do not need x-dir average
% 1 represents take 1 in x-direction

    u=reshape(uu(:,1),size_y,60)';
    v=reshape(vv(:,1),size_y,60)';
    w=reshape(ww(:,1),size_y,60)';

    umean(:,i) = mean(u,2);% average over y-direction
    vmean(:,i) = mean(v,2);
    wmean(:,i) = mean(w,2);
    
    umean_mid(i) = mean(umean(:,i));

    uflux(i)=sum(umean(:,i))*dz*dy;
    vflux(i)=sum(vmean(:,i))*dz*dx;
    wflux(i)=sum(wmean(:,i))*dx*dy;
    
   % u',v',w' are defined as u - <u>, where <.> is average over x,y direction 
    udiff(:,:,i) = u - umean(:,i)*ones(1,size_y);
    vdiff(:,:,i) = v - vmean(:,i)*ones(1,size_y);
    wdiff(:,:,i) = w - wmean(:,i)*ones(1,size_y);
    
   

end

%     figure(1)
%     clf
%     subplot(3,1,1)
%     plot(y,u_surf(:,1)/umean_mid(i),'LineWidth',2)
%     hold on
%     plot(y,u_mid(:,1)/umean_mid(i),'--','LineWidth',2)
%     plot(y,u_bot(:,1)/umean_mid(i),'.','LineWidth',2)
%     axis([min(y) max(y) 0.7 1.4])
%     ylabel('U/U_{mean}','FontSize',10);
%     %legend('Near surface','Middle depth','Near bottom')
%     title(['U (Time = ',t_str(i,:),' s)'])
%      
%     subplot(3,1,2)
%     plot(y,v_surf(:,1)/umean_mid(i),'LineWidth',2)
%     hold on
%     plot(y,v_mid(:,1)/umean_mid(i),'--','LineWidth',2)
%     plot(y,v_bot(:,1)/umean_mid(i),'.','LineWidth',2)
%     axis([min(y) max(y) -0.15 0.15])
%     ylabel('V/U_{mean}','FontSize',10);
%     %legend('Near surface','Middle depth','Near bottom')
%     title(['V (Time = ',t_str(i,:),' s)']) 
%     
%     
%     subplot(3,1,3)
%     plot(y,w_surf(:,1)/umean_mid(i),'LineWidth',2)
%     hold on
%     plot(y,w_mid(:,1)/umean_mid(i),'--','LineWidth',2)
%     plot(y,w_bot(:,1)/umean_mid(i),'.','LineWidth',2)
%     axis([min(y) max(y) -0.15 0.05])
%     ylabel('W/U_{mean}','FontSize',10);
%     %legend('Near surface','Middle depth','Near bottom')
%     title(['W (Time = ',t_str(i,:),' s)']) 
%     xlabel('Y (m)','FontSize',10);


udiff_max = max(max(max(udiff)));
udiff_min = min(min(min(udiff)));
    
vdiff_max = max(max(max(vdiff)));
vdiff_min = min(min(min(vdiff)));
    
wdiff_max = max(max(max(wdiff)));
wdiff_min = min(min(min(wdiff)));

umn_max = max(max(umean));
umn_min = min(min(umean));

vmn_max = max(max(vmean));
vmn_min = min(min(vmean));

wmn_max = max(max(wmean));
wmn_min = min(min(wmean));

%return  

mov1(1:nFrames) = struct('cdata', [],'colormap', []);
mov2(1:nFrames) = struct('cdata', [],'colormap', []);
mov3(1:nFrames) = struct('cdata', [],'colormap', []);

%% plot 

for i=1:nFrames;
%for i=7;   
    %% velocity difference contour
    figure(1);
    clf;
    subplot(3,1,1)
    contourf(y,z,udiff(:,:,i)/umean_mid(i),10);
   % caxis([udiff_min udiff_max]/umean_mid(i))
    caxis([-0.35 0.35])
    colorbar
    ylabel('Z (m)','FontSize',10);
    title(['U (Time = ',t_str(i,:),' s)'])
    
    subplot(3,1,2)
    contourf(y,z,vdiff(:,:,i)/umean_mid(i),10);
    %caxis([udiff_min udiff_max]/umean_mid(i))
    caxis([-0.35 0.35])
    colorbar
    ylabel('Z (m)','FontSize',10);
    title(['V (Time = ',t_str(i,:),' s)'])

    subplot(3,1,3)
    contourf(y,z,wdiff(:,:,i)/umean_mid(i),10);
    %caxis([udiff_min udiff_max]/umean_mid(i))
    caxis([-0.35 0.35])
    colorbar
    xlabel('Y (m)','FontSize',10)
    ylabel('Z (m)','FontSize',10);
    title(['W (Time = ',t_str(i,:),' s)'])
    
    print -djpeg -r300 uvw_contour2

    set(gca,'nextplot','replacechildren');
    mov(i) = getframe(gcf);
    
   % mean velocity plot

   
    figure(2)
    clf
    subplot(3,1,1)
    plot(umean(:,i)/umean_mid(i),z,'LineWidth',2)
    axis([0.5 2.5 -len_z -dz])
    ylabel('Z (m)','FontSize',10);
    title(['U (Time = ',t_str(i,:),' s)'])
    
    subplot(3,1,2)
    plot(vmean(:,i)/umean_mid(i),z,'LineWidth',2)
   % axis([vmn_min vmn_max -len_z -dz])
    axis([-1e-1 1e-1 -len_z -dz])
    ylabel('Z (m)','FontSize',10);
    title(['V (Time = ',t_str(i,:),' s)']) 
    
    subplot(3,1,3)
    plot(wmean(:,i)/umean_mid(i),z,'LineWidth',2)
   % axis([wmn_min wmn_max -len_z -dz])
    axis([-1e-3 1e-3 -len_z -dz])
    xlabel('Y (m)','FontSize',10)
    ylabel('Z (m)','FontSize',10);
    title(['W (Time = ',t_str(i,:),' s)'])
    
    set(gca,'nextplot','replacechildren');
    mov2(i) = getframe(gcf);
    %print -djpeg -r300 uvw_mean_profile2
    % surface elevation plot 
    
    figure(3)
    clf
    plot(y,eta(:,1),'LineWidth',2)
    hold on
    plot(y,y*0,'k--','LineWidth',2)
    axis([dy len_y -0.05 0.05])
    xlabel('Y (m)','FontSize',10);
    ylabel('Eta (m)','FontSize',10);

    set(gca,'nextplot','replacechildren');
    mov3(i) = getframe(gcf);

end
%return
movie2avi(mov, 'uvw1hr.avi', 'FPS',2,'compression', 'Cinepak','quality',50,'keyframe',1);
movie2avi(mov2, 'uvwmean1hr.avi', 'FPS',2,'compression','Cinepak','quality',50,'keyframe',1);
movie2avi(mov3, 'eta1hr.avi', 'FPS',2,'compression', 'Cinepak','quality',50,'keyframe',1);

% plot mean U profile comparison over time 

% figure(4);
% 
% u0=load('uprofstdy.mat'); % should be 50 hour run
% u0=u0.u;
% plot(u0,z,'LineWidth',2)
% hold on
% plot(umean(:,20),z,'r--','LineWidth',2)
% plot(umean(:,40),z,'k--','LineWidth',2)
% plot(umean(:,60),z,'g--','LineWidth',2)
% xlabel('U (m/s)','FontSize',14)
% ylabel('Z (m)','FontSize',14)
% h_legend = legend('No LC','LC 20 min','LC 40 min','LC 60 min');
% set(h_legend,'FontSize',14)


% plot flux in U, V, W direction

figure(5)
subplot(3,1,1)
plot(time,uflux./umean_mid/dy/len_z,'LineWidth',2)
ylabel('x-Flux','FontSize',10)

subplot(3,1,2)
plot(time,vflux./umean_mid/dx/len_z,'LineWidth',2)
axis([min(time) 4000 -2e-5 2e-5])
ylabel('y-Flux ','FontSize',10)

subplot(3,1,3)
plot(time,wflux./umean_mid/dy/dx/size_z,'LineWidth',2)
axis([min(time) 4000 -2e-5 2e-5])
xlabel('Time (s)','FontSize',10)
ylabel('z-Flux ','FontSize',10)















