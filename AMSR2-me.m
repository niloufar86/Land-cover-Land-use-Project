%%% Reading AMSR-2 0.25 deg gridded Tb product in H5 format      %%%%%%%%%%
clc;
close all;
clear all;
%
[ilat cellcntr icells box flat flon dlont thismax iind jind]=textread('p25ancil.out',...
    '%d%d%d%d%f%f%f%f%d%d');
P2=load('p25ancil2.out');
cellN = load('LandcellN.dat');
%
dlontb=P2(:,3);
ncells=P2(:,1);
%

lat1 = 0.125:0.25:179.875;
lat1=lat1';
lat2=repmat(lat1,1,1440);
lat2=180-lat2;
lon1 = 0.125:0.25:359.875;
lon2=repmat(lon1,720,1);
%
ff = [06; 07; 10; 18; 23; 36; 89];

%%%%%% READING AMSR-2 HDF FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for f = 1:7;
        freq = ff(f);
        
dirName=['/Users/student/Documents/NASA Project/AMSR2/2014.01/TB' num2str(ff(f), '%02d') 'GHz'];
files=dir(fullfile(dirName,'GW1AM2*_01D_EQMA*.h5'));
file_list={files.name}';
nfile=numel(file_list);

%
for day=1:nfile;
Fname = fullfile(dirName,file_list{day});
s = hdf5info(Fname);
tbh = double(hdf5read(s.GroupHierarchy.Datasets(1)));
tbv = double(hdf5read(s.GroupHierarchy.Datasets(2)));
%
tbh = tbh'*0.01;
tbv = tbv'*0.01;
tbh(tbh > 400 | tbh < 100) = NaN;
tbv(tbv > 400 | tbv < 100) = NaN;
%%%% Reprojecting into a matrix T of size 660066  %%%%%%%%%%%%%%
%
TTA1=zeros(660066,1);
CTA1=zeros(660066,1);
T1=zeros(660066,1);
TTA2=zeros(660066,1);
CTA2=zeros(660066,1);
T2=zeros(660066,1);
for i=1: 720*1440;
                 ieqlat=floor((lat2(i)/0.25)+1);
                 ieqlon=floor((lon2(i)/dlontb(ieqlat)) +1);
                 ibox=ncells(ieqlat)+ieqlon;
                 if (isfinite(tbh(i))>0)
                 TTA1(ibox)=tbh(i)+TTA1(ibox);
                 CTA1(ibox)=CTA1(ibox)+1;
                 end;
                 if (isfinite(tbv(i))>0)
                 TTA2(ibox)=tbv(i)+TTA2(ibox);
                 CTA2(ibox)=CTA2(ibox)+1;
                 end;
        end;
       T1= TTA1./CTA1;   
       T2= TTA2./CTA2; 
%%%% Extracting Land only data (Size should be 177499) %%%%%%%%%%%%%%%%%%      
      TB_H(:,day)=T1(cellN);
      TB_V(:,day)=T2(cellN);
end;
    
save(['/Users/student/Documents/NASA Project/' num2str(ff(f), '%02d') 'H_ASC_Jan2014.mat'],'TB_H');
save(['/Users/student/Documents/NASA Project/' num2str(ff(f), '%02d') 'V_ASC_Jan2014.mat'],'TB_V');

 end;     

%%%%%%%%%%%  Restructing for plotting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TT=zeros(660066,1);
TT(cellN)=TB_H(:,1);
TT(TT == 0) = NaN;
% %
 mtx=zeros(1440,720);
for i=1:1440*720
    mtx(i)= TT(box(i));
end;
mtx1=mtx;
mtx1(1:720,:)=mtx(721:1440,:);
mtx1(721:1440,:)=mtx(1:720,:);
% %%%%%%%%%%%%%%%% PLOTTING FIGURE %%%%%%%%%%%%%%%%%%%%%%%%
figure;
imagesc(flipud(mtx'));
%caxis([200 330]);
colormap(jet);
colorbar;
% %%%%%%%%%%%%%%%%%%%%%%% THE END %%%%%%%%%%%%%%%%%%%%%%%%