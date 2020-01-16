clc;
close all;
clear all;
%
[ilat cellcntr icells box flat flon dlont thismax iind jind]=textread(...
    'p25ancil.out','%d%d%d%d%f%f%f%f%d%d');
P2=load('p25ancil2.out');
cellN = load('LandcellN.dat');
%
dlontb=P2(:,3);
ncells=P2(:,1);

%%%%%% READING MODIS HDF FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirName=['/Users/student/Documents/NASA Project/MODIS/2014.01'];
files=dir(fullfile(dirName,'MYD11C1*.hdf'));
file_list={files.name};
nfile=numel(file_list);


%%%  Lat & Lon should be positive values only   %%%%%%%%%%
%%%  Matrix size of Lat & Lon should be same as LST size %%%%%%%%%%
lat1 = 0.025:0.05:179.975;
lat1=lat1';
lat2=repmat(lat1,1,7200);
lat2=180-lat2;
lon1 = 0.025:0.05:359.975;
lon2=repmat(lon1,3600,1);
%
for day=1:nfile
Fname = fullfile(dirName,file_list{day});
data1 = hdfread(Fname, '/MODIS_CMG_3MIN_LST/Data Fields/LST_Day_CMG', 'Index', {[1  1],[1  1],[3600  7200]});
lst = double(data1);
lst = lst *0.02;
lst(lst > 655.0 | lst < 75.0) = NaN;
%
%%% Converting -180 to 180 Longitude into 0 to 360 format %%%%%%%%%
lst1=lst;
lst1(:,1:3600)=lst(:,3601:7200);
lst1(:,3601:7200)=lst(:,1:3600);
%
%%%% Reprojecting into 0.25 deg grid resolution %%%%%%%%%%%%%%%%%%%%
%%%%     Size of the final matrix T should be 660066  %%%%%%%%%%%%%%
%
TTA=zeros(660066,1);
CTA=zeros(660066,1);
for i=1: 3600*7200
                 ieqlat=floor((lat2(i)/0.25)+1);
                 ieqlon=floor((lon2(i)/dlontb(ieqlat)) +1);
                 ibox=ncells(ieqlat)+ieqlon;
                 if (isfinite(lst1(i))>0)
                 TTA(ibox)=lst1(i)+TTA(ibox);
                 CTA(ibox)=CTA(ibox)+1;
                 end;
        end;
       T= TTA./CTA;   
%%%% Extracting Land only data (Size should be 177499) %%%%%%%%%%%%%%%%%%      
      LANDT(:,day)=T(cellN);
end;
save(['/Users/student/Documents/NASA Project/LST_Jan2014.mat'],'LANDT');

%%%%%%%%%%%  Restructing for plotting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TT=zeros(660066,1);
TT(cellN)=LANDT(:,1);
TT(TT == 0) = NaN;
% 
mtx=zeros(1440,720);
for i=1:1440*720
mtx(i)= TT(box(i));
end;
% %%%%%%%%%%%%%%% PLOTTING FIGURE %%%%%%%%%%%%%%%%%%%%%%%%
figure;
load coast;
maplegend=[1/.25 90 360];
ah=axesm ('MapProjection','eqdcylin',...
'MapLatLimit',[-90 90],'MapLonLimit',[-180 180],'Frame','on','Grid','on');
%cmap=colormap(jet(80));
%colormap(cmap)
colormap(jet(80))
clear mlat mlon
[mlat,mlon]=meshgrat(mtx',maplegend);
surfm(mlat,mlon,mtx','Facecolor','interp');
geoshow(lat, long);
contourcmap([240:10:320], 'jet', 'colorbar', 'on','Location' ,'horizontal');