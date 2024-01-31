function [simImg,mask,parMap,smap,S0,ID,aif,T10] = gen_DRO(data,option)


idx = round(length(data)*rand());
% idx = 12;

par_var = 0.1;
par_var_t = 0.2;
t = linspace(0,150,22);

mask_ = data{idx}.mask;
if ~isfield(option,'AIF')
    aif = data{idx}.AIF;
else
    aif = option.AIF;
end
S0 = data{idx}.S0;
% ID = data{idx}.ID;
smap = double(data{idx}.smap);
% disp(ID);

if ~isfield(option,'B1')
    B1 = ones(size(S0));
end
if isfield(option,'parMap')
    parMap = option.parMap;
end

% 1-gland, 2-benign, 3-malig, 4-muscle, 5-skin, 6-liver, 7-heart, 8-vasc
mask.glandular = mask_==1;
mask.benign = mask_==2;
mask.malignant = mask_==3;
mask.muscle = mask_==4;
mask.skin = mask_==5;
mask.liver = mask_==6;
mask.heart = mask_==7;
mask.vascular = mask_==8;

% 
% if ~isfield(mask,'heart')
%     mask.heart = logical(zeros(size(S0,1),size(S0,2)));
% end
% if ~isfield(mask,'liver')
%     mask.liver = logical(zeros(size(S0,1),size(S0,2)));
% end
% if ~isfield(mask,'skin')
%     mask.skin = logical(zeros(size(S0,1),size(S0,2)));
% end
% if ~isfield(mask,'muscle')
%     mask.muscle = logical(zeros(size(S0,1),size(S0,2)));
% end
% if ~isfield(mask,'benign')
%     mask.benign = logical(zeros(size(S0,1),size(S0,2)));
% end
% if ~isfield(mask,'malignant')
%     mask.malignant = logical(zeros(size(S0,1),size(S0,2)));
% end
% if ~isfield(mask,'heart_blood')
%     mask.heart_blood = mask.heart;
% end


nbase = find((aif<0.15)&(t<300));
nbase = nbase(end);


t_1s = 0:t(end);
aifci = interp1(t,aif,t_1s,'pchip');

[nx,ny] = size(S0);
%logIdx = gen_downSample_logIdx(t,ti);

mask_inner = (mask.heart | mask.liver);
se = strel('disk',4);
mask_inner = imdilate(mask_inner,se);

T1.glandular = 1.324;
T1.malignant = 1.5;
T1.benign = 1.4;
T1.liver = 0.81;
% T1.heart = 1.68;
T1.heart = 0.81;
T1.muscle = 1.41;
T1.skin = 0.85; %Check
% T1.skin = 1;
T1.vasc = 1.93;

T10 = zeros(nx,ny);
T10(mask.glandular) = T1.glandular;
T10(mask.malignant) = T1.malignant;
T10(mask.benign) = T1.benign;
T10(mask.liver) = T1.liver;
T10(mask.heart) = T1.heart;
T10(mask.muscle) = T1.muscle;
T10(mask.skin) = T1.skin;
T10(mask.vascular) = T1.heart;
% T10(T10==0) = 1e-8;
T10(T10==0) = 1;
temp = T10;
% temp(~mask_inner) = 1;
temp = imgaussfilt(temp,10);
T10(mask_inner) = temp(mask_inner);


if ~exist('parMap')
p0.glandular = [0,0; 0.010,0.077; 0.1/60,0.115/60; 0.01/60,0.043/60];
p0.malignant = [0.101,0.3;0.131,0.256;0.259/60,1.032/60;0.0434/60,1.98/60];
p0.benign = [0.141,0.3;0.011,0.190;0.116/60,0.228/60;0.05/60,1.056/60];
p0.muscle = [0,0; 0.010,0.101; 0.1/60,0.118/60; 0.011/60,0.069/60];
p0.skin = [0,0; 0.039,0.125; 0.1/60,0.151/60; 0.01/60,0.019/60]; 

p0.liver = [0.1,0.5;0.353,0.5;0.433/60,1.227/60;1.990/60,2/60];
p0.heart = [0.148,0.300;0.214,0.373;2/60,2/60;0.404/60,1.224/60];
p0.vascular = [0,0;0.3,0.3;2/60,2/60;0/60,0/60];

gland_ktrans = [0.01,0.0352]./60;
malig_ktrans = [0.0412,0.385]./60;
benign_ktrans = [0.0453,0.143]./60;
muscle_ktrans = [0.011,0.05]./60;
skin_ktrans = [0.009,0.017]./60;
liver_ktrans = [0.412,0.979]./60;
heart_ktrans = [0.365,0.810]./60;


fun_ktrans = @(x) x(3).*(1-exp(-x(4)./x(3)));
p0_glandular = p0.glandular(:,1)+(p0.glandular(:,2)-p0.glandular(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_glandular);
while ((par_ktrans)<gland_ktrans(1)||(par_ktrans)>gland_ktrans(2))
    p0_glandular = p0.glandular(:,1)+(p0.glandular(:,2)-p0.glandular(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_glandular);
end
p0_malignant = p0.malignant(:,1)+(p0.malignant(:,2)-p0.malignant(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_malignant);
while ((par_ktrans)<malig_ktrans(1)||(par_ktrans)>malig_ktrans(2))
    p0_malignant = p0.malignant(:,1)+(p0.malignant(:,2)-p0.malignant(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_malignant);
end
p0_benign = p0.benign(:,1)+(p0.benign(:,2)-p0.benign(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_benign);
while ((par_ktrans)<benign_ktrans(1)||(par_ktrans)>benign_ktrans(2))
    p0_benign = p0.benign(:,1)+(p0.benign(:,2)-p0.benign(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_benign);
end
p0_muscle = p0.muscle(:,1)+(p0.muscle(:,2)-p0.muscle(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_muscle);
while ((par_ktrans)<muscle_ktrans(1)||(par_ktrans)>muscle_ktrans(2))
    p0_muscle = p0.muscle(:,1)+(p0.muscle(:,2)-p0.muscle(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_muscle);
end
p0_skin = p0.skin(:,1)+(p0.skin(:,2)-p0.skin(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_skin);
while ((par_ktrans)<skin_ktrans(1)||(par_ktrans)>skin_ktrans(2))
    p0_skin = p0.skin(:,1)+(p0.skin(:,2)-p0.skin(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_skin);
end
p0_heart = p0.heart(:,1)+(p0.heart(:,2)-p0.heart(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_heart);
while ((par_ktrans)<heart_ktrans(1)||(par_ktrans)>heart_ktrans(2))
    p0_heart = p0.heart(:,1)+(p0.heart(:,2)-p0.heart(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_heart);
end
p0_liver = p0.liver(:,1)+(p0.liver(:,2)-p0.liver(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_liver);
while ((par_ktrans)<liver_ktrans(1)||(par_ktrans)>liver_ktrans(2))
    p0_liver = p0.liver(:,1)+(p0.liver(:,2)-p0.liver(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_liver);
end

p0_vascular = p0.vascular(:,1)+(p0.vascular(:,2)-p0.vascular(:,1)).*rand(4,1);

parMap = zeros(nx,ny,4);

for i = 1:4
    temp = zeros(nx,ny);
    randMap = p0_glandular(i)*((1-par_var)+(par_var*2)*rand(nx,ny));
    temp(mask.glandular) = randMap(mask.glandular);
    
    randMap = p0_malignant(i)*((1-par_var_t)+(par_var_t*2)*rand(nx,ny));
    temp(mask.malignant) = randMap(mask.malignant);
    
    randMap = p0_benign(i)*((1-par_var_t)+(par_var_t*2)*rand(nx,ny));
    temp(mask.benign) = randMap(mask.benign);
    
    %temp(mask.fat) = p0.fat(i);
    randMap = p0_liver(i)*((1-par_var)+(par_var*2)*rand(nx,ny));
    temp(mask.liver) = randMap(mask.liver);
    
    %temp(mask.myocard) = p0.myocard(i);
    randMap = p0_muscle(i)*((1-par_var)+(par_var*2)*rand(nx,ny));
    temp(mask.muscle) = randMap(mask.muscle);
    
    randMap = p0_skin(i)*((1-par_var)+(par_var*2)*rand(nx,ny));
    temp(mask.skin) = randMap(mask.skin);
    %temp = temp.*(0.8+(0.4*rand(nx,ny)));
    
    randMap = p0_vascular(i)*ones(nx,ny);
    temp(mask.vascular) = randMap(mask.vascular);
    

    randMap = p0_heart(i)*((1-par_var)+(par_var*2)*rand(nx,ny));
    temp(mask.heart) = randMap(mask.heart);
    parMap(:,:,i) = temp;
end
end

aifci_1s = zeros(nx,ny,length(t_1s));
[rIdx,cIdx] = find(mask.liver==1); %10s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci([1:10,1:end-10]);
end

[rIdx,cIdx] = find(mask.glandular==1);%15s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci([1:15,1:end-15]);
end

[rIdx,cIdx] = find((mask.vascular|mask.heart)==1); %0s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci;
end


[rIdx,cIdx] = find(mask.malignant==1); %3s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci([1:3,1:end-3]);
end
[rIdx,cIdx] = find(mask.benign==1); %8s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci([1:8,1:end-8]);
end

[rIdx,cIdx] = find(mask.muscle==1); %7s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci([1:7,1:end-7]);
end
[rIdx,cIdx] = find(mask.skin==1); %10s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci([1:5,1:end-5]);
end
% bgd_mask = parMap(:,:,4)==0;


ti = 0:0.1:t(end);

logIdx = zeros(1,length(ti));
start_idx = 1;
for i = 1:length(t)
    for j = start_idx:length(ti)
        if t(i)<=ti(j)
            logIdx(j) = 1;
            start_idx = j;
            break
        end
    end
end
logIdx = logical(logIdx);

aifci_Map = zeros(nx,ny,length(ti));
[rIdx,cIdx] = find(parMap(:,:,2)>0);
for i = 1:length(rIdx)
    temp_aif = squeeze(aifci_1s(rIdx(i),cIdx(i),:));
    aifci_Map(rIdx(i),cIdx(i),:) = interp1(t_1s,temp_aif,ti,'pchip');
end


parMap(parMap==0) = 1e-8;
% mask_inner = edge(mask.heart | mask.liver);

for i = 1:4
%     parMap(:,:,i) = imgaussfilt(parMap(:,:,i));
    par_temp = imgaussfilt(parMap(:,:,i),1);
    temp = parMap(:,:,i);
    temp(~mask_inner) = 0;
    temp = imgaussfilt(temp,20);
    % kernel = ones(windowSize)/windowSize^2;
    % temp = imfilter(temp,kernel);
    par_temp(mask_inner) = temp(mask_inner);
    parMap(:,:,i) = par_temp;
end


ve = parMap(:,:,1);
vp = parMap(:,:,2);
fp = parMap(:,:,3);
ktrans = parMap(:,:,4);


Ce = zeros(size(parMap,1),size(parMap,2),size(aifci_Map,3));
Cp = zeros(size(parMap,1),size(parMap,2),size(aifci_Map,3));

for i = 2:size(aifci_Map,3)
    dt = ti(i)-ti(i-1);
    dcp = fp.*aifci_Map(:,:,i-1) - (fp+ktrans).*Cp(:,:,i-1) + ktrans.*Ce(:,:,i-1);
    dce = ktrans.*Cp(:,:,i-1)-ktrans.*Ce(:,:,i-1);
    Cp(:,:,i) = Cp(:,:,i-1) + dcp.*dt./vp;
    Ce(:,:,i) = Ce(:,:,i-1) + dce.*dt./ve;
end

cts = Cp.*repmat(vp,[1,1,length(ti)]) + Ce.*repmat(ve,[1,1,length(ti)]);
cts = cts(:,:,logIdx);

% nan_mask = sum(isnan(cts),3);
% [rIdx,cIdx] = find(nan_mask);
% disp(['Found ',num2str(length(rIdx)),' NaN voxels']);
% for i = 1:length(rIdx)
%     cts(rIdx(i),cIdx(i),:) = 0;
% end

cts_tcm = cts;


ve = 1;
vp = parMap(:,:,2);
fp = parMap(:,:,3);
ktrans = parMap(:,:,4);

Ce = zeros(size(parMap,1),size(parMap,2),size(aifci_Map,3));
Cp = zeros(size(parMap,1),size(parMap,2),size(aifci_Map,3));

for i = 2:size(aifci_Map,3)
    dt = ti(i)-ti(i-1);
    dcp = fp.*aifci_Map(:,:,i-1) - (fp+ktrans).*Cp(:,:,i-1);
    dce = ktrans.*Cp(:,:,i-1);
    Cp(:,:,i) = Cp(:,:,i-1) + dcp.*dt./vp;
    Ce(:,:,i) = Ce(:,:,i-1) + dce.*dt./ve;
end

cts = Cp.*repmat(vp,[1,1,length(ti)]) + Ce.*repmat(ve,[1,1,length(ti)]);
cts = cts(:,:,logIdx);

% nan_mask = sum(isnan(cts),3);
% [rIdx,cIdx] = find(nan_mask);
% disp(['Found ',num2str(length(rIdx)),' NaN voxels']);
% for i = 1:length(rIdx)
%     cts(rIdx(i),cIdx(i),:) = 0;
% end

cts_epm = cts;



% cts(mask_vasc) = aif_mat(mask_vasc);

ve = 1;
vp = parMap(:,:,2);
ktrans = parMap(:,:,4);

Ce = zeros(size(parMap,1),size(parMap,2),size(aifci_Map,3));
Cp = zeros(size(parMap,1),size(parMap,2),size(aifci_Map,3));

for i = 2:size(aifci_Map,3)
    dt = ti(i)-ti(i-1);
%     dcp = fp.*aifci_Map(:,:,i-1) - (fp+ktrans).*Cp(:,:,i-1);
    dce = ktrans.*aifci_Map(:,:,i-1);
%     Cp(:,:,i) = Cp(:,:,i-1) + dcp.*dt./vp;
    Ce(:,:,i) = Ce(:,:,i-1) + dce.*dt./ve;
end

cts = aifci_Map.*repmat(vp,[1,1,length(ti)]) + Ce;
cts = cts(:,:,logIdx);

% nan_mask = sum(isnan(cts),3);
% [rIdx,cIdx] = find(nan_mask);
% disp(['Found ',num2str(length(rIdx)),' NaN voxels']);
% for i = 1:length(rIdx)
%     cts(rIdx(i),cIdx(i),:) = 0;
% end

cts_eTofts = cts;


mask_epm = mask.glandular | mask.skin | mask.muscle;
mask_epm = repmat(mask_epm,[1,1,size(cts,3)]);
mask_eTofts = mask.vascular;
mask_eTofts = repmat(mask_eTofts,[1,1,size(cts,3)]);

cts_tcm(mask_epm) = cts_epm(mask_epm);
cts_tcm(mask_eTofts) = cts_eTofts(mask_eTofts);
cts = cts_tcm;

nan_mask = sum(isnan(cts),3);
[rIdx,cIdx] = find(nan_mask);
disp(['Found ',num2str(length(rIdx)),' NaN voxels']);


TR = 4.87e-3;
theta = 10 * pi /180;

theta = theta .* B1;

r1 = 4.3;
T1_t = 1./(1./repmat(T10,[1,1,size(cts,3)]) + r1*cts);

Eh = (1 - exp(-TR./T1_t)).*sin(theta)./ (1-(exp(-TR./T1_t).*cos(theta)));
Eh = Eh ./ repmat(mean(Eh(:,:,1:nbase),3),[1,1,size(Eh,3)]);
%Eh(mask_enh) = 1e-8;


simImg = Eh .* repmat(S0,[1,1,size(Eh,3)]);
% aif_map = aifci_Map(:,:,logIdx);

nan_mask = sum(isnan(simImg),3);
[rIdx,cIdx] = find(nan_mask);
disp(['Sig-Found ',num2str(length(rIdx)),' NaN voxels']);
if nnz(nan_mask)>0
    close all;
    figure;
    imshow(abs(S0),[]);colormap(gray);
    hold on;freezeColors();
    h = imshow(nan_mask>0,[]);colormap(jet);
    set(h,'AlphaData',nan_mask>0);

    keyboard;
end

%simImg = simImg(:,:,logIdx);
%fsimImg = simImg / max(simImg(:));




