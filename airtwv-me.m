clc;
clear all;
close all;
%
[ilat cellcntr icells box flat flon dlont thismax iind jind]=textread('p25ancil.out',...
    '%d%d%d%d%f%f%f%f%d%d');
P2=load('p25ancil2.out');
cellN = load('LandcellN.dat');
%
dlontb=P2(:,3);
ncells=P2(:,1);

%%% Note that Lat & Lon should be positive values only   %%%%%%%%%%
lat1 = 0.125:0.25:179.875;
lat1=lat1';
lat2=repmat(lat1,1,1440);
lat2=180-lat2;
lon1 = 0.125:0.25:359.875;
lon2=repmat(lon1,720,1);
%
%%%%%% READING AIRS HDF FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirName=['/Users/student/Documents/NASA Project/AIRS/2014.01'];
files=dir(fullfile(dirName,'AIRS.2014.01.*.hdf'));
file_list={files.name};
nfile=numel(file_list);
%

for day=1:nfile
FILE_NAME = fullfile(dirName,file_list{day});

data = hdfread(FILE_NAME, '/ascending/Data Fields/TotH2OVap_A', 'Index', {[1  1],[1  1],[180  360]});

%

data=double(data);

data(data == -9999.0) = NaN;

%
%%% Converting -180 to 180 Longitude into 0 to 360 format %%%%%%%%%
airtwv=data;
airtwv(:,1:180)=data(:,181:360);
airtwv(:,181:360)=data(:,1:180);
%
%%%% Resampling into 0.25 degree grid resolution %%%%%%%%%%%%%%%%%
airtwv_p25=imresize(airtwv,[720 1440],'nearest');
%
%%%% Reprojecting into a matrix T of size 660066  %%%%%%%%%%%%%%
%
TTA=zeros(660066,1);
CTA=zeros(660066,1);
TWV=zeros(660066,1);
for i=1: 720*1440
                 ieqlat=floor((lat2(i)/0.25)+1);
                 ieqlon=floor((lon2(i)/dlontb(ieqlat)) +1);
                 ibox=ncells(ieqlat)+ieqlon;
                 if (isfinite(airtwv_p25(i))>0)
                 TTA(ibox)=airtwv_p25(i)+TTA(ibox);
                 CTA(ibox)=CTA(ibox)+1;
                 end;
        end;
       TWV= TTA./CTA;   
%%%% Extracting Land only data (Size should be 177499) %%%%%%%%%%%%%%%%%%      
      LandTWV(:,day)=TWV(cellN);
end;


save(['/Users/student/Documents/NASA Project/TotH2OVap_A_Jan2014.mat'],'LandTWV');

%%%%%%%%%%%  Restructing for plotting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TT=zeros(660066,1);
TT(cellN)=LandTWV(:,1);
TT(TT == 0) = NaN;
% %
mtx=zeros(1440,720);
for i=1:1440*720
  mtx(i)= TT(box(i));
end;
% %%%%%%%%%%%%%%%% PLOTTING FIGURE %%%%%%%%%%%%%%%%%%%%%%%%
figure;
imagesc(flipud(mtx'));
colormap(jet);
colorbar;
%%%%%%%%%%%%%%%%%%%%%%% THE END %%%%%%%%%%%%%%%%%%%%%%%%