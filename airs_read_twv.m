%%%%% To extract Total integrated column Water Vapor (kg/m^2) from AIRS %%
%%%%%%%%%%%%%%%  Satya Prakash, NYCCT, October 30, 2016  %%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
%
[ilat cellcntr icells box flat flon dlont thismax iind jind]=textread(...
    '/Users/hamidnorouzi/Documents/MATLAB/MODIS/p25ancil.out','%d%d%d%d%f%f%f%f%d%d');
P2=load('/Users/hamidnorouzi/Documents/MATLAB/MODIS/p25ancil2.out');
cellN = load('/Users/hamidnorouzi/Documents/MATLAB/MODIS/LandcellN.dat');
%
dlontb=P2(:,3);
ncells=P2(:,1);
%
mm = [01; 02; 03; 04; 05; 06; 07; 08; 09; 10; 11; 12];
mmm = ['JAN'; 'FEB'; 'MAR'; 'APR'; 'MAY'; 'JUN'; 'JUL'; 'AUG'; 'SEP'; 'OCT'; 'NOV'; 'DEC'];
%%% Note that Lat & Lon should be positive values only   %%%%%%%%%%
lat1 = 0.125:0.25:179.875;
lat1=lat1';
lat2=repmat(lat1,1,1440);
lat2=180-lat2;
lon1 = 0.125:0.25:359.875;
lon2=repmat(lon1,720,1);
%
for m = 1:7
    mon = mmm(m,:)
    clear LANDT T TTA CTA
%%%%%% READING AIRS HDF FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirName=['/Volumes/G-RAIDT/AIRS/2016/2016.' num2str(mm(m), '%02d')];
files=dir(fullfile(dirName,'AIRS.2016.*'));
file_list={files.name}';
nfile=numel(file_list)
%
%%%% Open the HDF-EOS2 Grid File.
for day=1:nfile
FILE_NAME = fullfile(dirName,file_list{day});

data = hdfread(FILE_NAME, '/ascending/Data Fields/TotH2OVap_A', 'Index', {[1  1],[1  1],[180  360]});
% file_id = hdfgd('open', FILE_NAME, 'rdonly');
% GRID_NAME='ascending';
% grid_id = hdfgd('attach', file_id, GRID_NAME);
% DATAFIELD_NAME='TotH2OVap_A';
% [data, fail] = hdfgd('readfield', grid_id, DATAFIELD_NAME, [], [], []);
%
data=double(data);
%[fillvalue,status] = hdfgd('getfillvalue',grid_id, DATAFIELD_NAME);
data(data == -9999.0) = NaN;
% data = data';
% hdfgd('detach', grid_id);
%
%%% Converting -180 to 180 Longitude into 0 to 360 format %%%%%%%%%
airt=data;   
airt(:,1:180)=data(:,181:360);
airt(:,181:360)=data(:,1:180);
%
%%%% Resampling into 0.25 degree grid resolution %%%%%%%%%%%%%%%%%
airt_p25=imresize(airt,[720 1440],'nearest');
%
%%%% Reprojecting into a matrix T of size 660066  %%%%%%%%%%%%%%
%
TTA=zeros(660066,1);
CTA=zeros(660066,1);
T=zeros(660066,1);
for i=1: 720*1440
                 ieqlat=floor((lat2(i)/0.25)+1);
                 ieqlon=floor((lon2(i)/dlontb(ieqlat)) +1);
                 ibox=ncells(ieqlat)+ieqlon;
                 if (isfinite(airt_p25(i))>0)
                 TTA(ibox)=airt_p25(i)+TTA(ibox);
                 CTA(ibox)=CTA(ibox)+1;
                 end;
        end;
       T= TTA./CTA;   
%%%% Extracting Land only data (Size should be 177499) %%%%%%%%%%%%%%%%%%      
      LANDT(:,day)=T(cellN);
end;

save(['/Volumes/G-RAIDT/Emissivity-AMSR2/TWV/2016/AIRS_TWV_ASC_' mon '2016.mat'],'LANDT');
end;
%%%%%%%%%%%  Restructing for plotting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TT=zeros(660066,1);
% TT(cellN)=LANDT(:,1);
% TT(TT == 0) = NaN;
% %
% mtx=zeros(1440,720);
% for i=1:1440*720
%     mtx(i)= TT(box(i));
% end;
% %%%%%%%%%%%%%%%% PLOTTING FIGURE %%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% imagesc(flipud(mtx'));
% colormap(jet);
% colorbar;
%%%%%%%%%%%%%%%%%%%%%%% THE END %%%%%%%%%%%%%%%%%%%%%%%%