close all;
clear all;

D1=0.25;
D2=0.05;
m1=719;
n1=1439;

maxm2=(m1*D1)/D2;
minm2=((m1*D1)-D1+D2)/D2;

maxn2=(n1*D1)/D2;
minn2=((n1*D1)-D1+D2)/D2;


%%%%% MODIS 0.05 degree 2014 LST data %%%%%


dirName=['/Users/student/Documents/NASA Project/MODIS/2014.01']
files=dir(fullfile(dirName,'MYD11C1*.hdf'));
file_list={files.name};
nfile=numel(file_list);
T=nan(3600,7200);

for day=1:nfile
Fname = fullfile(dirName,file_list{day});
dat = hdfread(Fname, '/MODIS_CMG_3MIN_LST/Data Fields/LST_Day_CMG', 'Index', {[1  1],[1  1],[3600  7200]});
lst = double(dat);
lst = lst *0.02;
lst(lst > 655.0 | lst < 75.0) = NaN;
T=cat(3,T,lst);
end;
T(:,:,1)=[];


A=T(minm2:maxm2,minn2:maxn2,:);
data1=reshape(A,[25,nfile]);
% a=find(isnan(sum(data1))==1);
% data1(:,a)=[];
data=data1';
figure;
imagesc(A(:,:,1))

[coeff1,score1,latent,tsquared,explained,mu1]=pca(data, 'algorithm','als');
LST=score1*coeff1'+repmat(mu1,31,1);
LST=LST';
LST1=reshape(LST',5,5,31)
LST1=real(LST1);
figure;
imagesc(LST1(:,:,1))

figure;
imagesc(LST1(:,:,1)-A(:,:,1));

figure;
imagesc(A(:,:,1));

%save(['/Users/student/Documents/NASA Project/PCA/PCA_MODIS_' mon '2014.mat'],'data','coeff1','score1','latent','tsquared','explained','mu1');

