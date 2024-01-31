function [kspace,traj,test] = gen_kspace_data(img,smap,spokes_per_frame,n_lvl)

img(img==0) = 1e-8;

addpath('nufft_toolbox/');
nx = size(img,1)*2;
ntviews = spokes_per_frame * size(img,3);
traj = Trajectory_GoldenAngle_GROG(ntviews,nx)/nx;
dcf = repmat(abs(linspace(-1,1,nx))',[1,ntviews]);
for t = 1:size(img,3)
    dcfu(:,:,t) = dcf(:,(t-1)*spokes_per_frame+1:t*spokes_per_frame);
    traju(:,:,t) = traj(:,(t-1)*spokes_per_frame+1:t*spokes_per_frame);
end

smap=double(smap/max(abs(smap(:))));
%     smap = ones(size(smap));
param.E = MCNUFFT(traju,dcfu,smap);

kspace = param.E * double(img);

kspace = kspace + n_lvl * randn(size(kspace));

kspace = kspace ./ repmat(reshape(sqrt(dcfu),[nx,spokes_per_frame,1,size(img,3)]),[1,1,size(smap,3),1]);

test = param.E' * kspace;
%     test = test./max(test(:));
kspace = permute(kspace,[1,2,4,3]); % Swapping time and coil dim
kspace = reshape(kspace,[nx,size(kspace,2)*size(kspace,3),size(kspace,4)]);

 





