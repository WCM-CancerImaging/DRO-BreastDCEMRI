
% Data collection

clc;clear;
addpath('./data');
for c = 1:53
    load(['sub',num2str(c),'.mat']);
    data{c}.AIF = aif;
    data{c}.mask = mask;
    data{c}.S0 = S0;
    data{c}.smap = smap;
end

%% Generate DRO

option = [];
[simImg,mask,parMap,smap,S0] = gen_DRO(data,option);

close all;
figure(100)
for i = 1:size(simImg,3)
    imshow(simImg(:,:,i),[0,max(simImg(:))],'InitialMagnification',1600);
    frame = getframe(100);
    img = frame2im(frame);
    [imind cm] = rgb2ind(img,256);
end

%% Display masks

disp_mask(mask,S0,1);

%% Radial k-space data generation

spokes_per_frame = 13;
n_lvl = 0.05;
[kspace,traj] = gen_kspace_data(simImg,smap,spokes_per_frame,n_lvl);


%% Radial reconstruction using BART
% NEED BART INSTALLATION TO RUN!

nt = floor(size(kspace,2)/spokes_per_frame);
kspace_trim = kspace(:,1:spokes_per_frame*nt,:);
traj_trim = traj(:,1:spokes_per_frame*nt);

[nx,ntview,ncoil] = size(kspace_trim);
kspace_dim = reshape(kspace_trim,[1,nx,spokes_per_frame,nt,ncoil]);
kspace_dim = permute(kspace_dim,[1,2,3,5,4]);
kspace_dim = reshape(kspace_dim,[1,nx,spokes_per_frame,ncoil,1,1,1,1,1,1,nt]);

clear traj_dim
traj_dim(1,:,:) = real(traj_trim);
traj_dim(2,:,:) = imag(traj_trim);
traj_dim = cat(1,traj_dim,zeros(1,nx,ntview));
traj_dim = traj_dim * (nx/2);
traj_dim = reshape(traj_dim,[3,nx,spokes_per_frame,ones(1,7),nt]);

smap_dim = reshape(smap,[size(smap,1),size(smap,2),1,ncoil]);

% TV Regularization, lamb=0.01
reco = bart('pics -S -RT:1024:0:0.01 -i100 -t',traj_dim,kspace_dim,smap_dim);

grasp_bart = (squeeze(abs(reco)));


