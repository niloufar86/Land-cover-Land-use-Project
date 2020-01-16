%%%% Computation of emissivity at AMSR-2 frequencies after computing 
%%%% atmospheric extinction, upwelling and downwelling Tbs from
%%%%   AIRS air temperaure (K) profile and TWV (kg/m^2) 
%%%%                      NYCCT: October 30, 2016
%%%%        NOTE:  DAY--> ASCENDING   &  NIGHT-->DESCENDING
%%%%               Modified on November 04, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
%
FF = [6.925; 7.30; 10.65; 18.7; 23.8; 36.5; 89.0];
FF1 = [06; 07; 10; 18; 23; 36; 89];
mm = [01; 02; 03; 04; 05; 06; 07; 08; 09; 10; 11; 12];
dd = [31; 29; 31; 30; 31; 30; 31; 31; 30; 31; 30; 31];
mmm = ['JAN'; 'FEB'; 'MAR'; 'APR'; 'MAY'; 'JUN'; 'JUL'; 'AUG'; 'SEP'; 'OCT'; 'NOV'; 'DEC']; 
UMU = cosd(55);  %AMSR-2 incidence angle
%
for id = 1:12
    clear mon no_day TWV AIRT SKINT WV Temp_Profile emis_h emis_v
    mon = mmm(id,:)
    no_day = dd(id);
load(['/Volumes/G-RAIDT/Emissivity-AMSR2/TWV/2016/AIRS_TWV_ASC_' mon '2016.mat']);
TWV = LANDT./10;  %%Conveting kg/m^2 to cm
%
load(['/Volumes/G-RAIDT/Emissivity-AMSR2/TEMP/2016/AIRS_TEMP_ASC_' mon '2016.mat']);
AIRT = LANDT;
%
load(['/Volumes/G-RAIDT/Emissivity-AMSR2/LST/2016/MODIS_LST_DAY_' mon '2016.mat']);
SKINT = LANDT;
%
for day=1:no_day
    display(day)
    
WV = TWV(:,day);
        
Temp_Profile = AIRT(:,:,day);
 
for f=1:7
    freq=FF(f)
    clear BTEMPV BTEMPH TBV TBH EMISV EMISH

    a = 0+(2.18*(10^-3))*freq+(3.9*10^-4)*freq^2;
    b = 9.73+(-8.92*(10^-2))*freq+(1.73*10^-4)*freq^2;

        Tlayers=[]; alfalayers=[]; alfatau=[]; alfa=[]; alfatauupw=[];
        alfatau=[]; Tauupw21=[]; Taudw2=[]; taudw=[]; Taudw21=[];
        Tauupw2=[]; tauupw=[]; alfa2=[]; T2=[]; T=[]; Tau=[]; Tup=[]; 
        Tdown=[]; no_temp=[]; alfataudw=[]; Teta=[]; w=[];

        for tt=1:12

        temp = squeeze(Temp_Profile(:,tt));
        no_temp = find(isnan(temp));
        T = temp; 
        T(no_temp) = 220;
        Teta = 300./T; % T must be in Kelvin

        w = (0.002166 *100.*WV)./T; %Converting TWV (cm) to concentration

        alfa = w.*(a.*(Teta.^(b)));

        alfa(no_temp) = 0;
        alfa(isnan(WV)) = 0;

        alfalayers = cat(3,alfalayers,alfa);
        alfatau = shiftdim(alfalayers,2);

        Tlayers = cat(3,Tlayers,T);

        end
c=1;

a2=size(alfalayers);

        for c = 1:a2(3) 
            
      if c==1
        taudw(:,:,c) = alfatau(c,:,:);
        tauupw(:,:,c) = trapz(alfatau);
      end

      if c>1 && c<a2(3)
        alfataudw = alfatau(1:c,:,:);
        alfatauupw = alfatau(c:a2(3),:,:);

        taudw(:,:,c) = trapz(alfataudw); 
        tauupw(:,:,c) = trapz(alfatauupw); 
      end

      if c==a2(3)
        taudw(:,:,c) = trapz(alfatau);
        tauupw(:,:,c) = alfatau(c,:,:);
      end

        end

%%%%********************** Integration ************************************
          alfa2 = shiftdim(alfalayers,2);
          T2 = shiftdim(Tlayers,2);
          
          Taudw2 = shiftdim(taudw,2);
          Tauupw2 = shiftdim(tauupw,2);
          Taudw21 = zeros(12,177499);
          Taudw21(:) = Taudw2(:);
          Tauupw21 = zeros(12,177499);
          Tauupw21(:) = Taudw2(:);
          
          Tdown = squeeze(trapz(T2.*alfa2.*exp(-Taudw21./UMU),1));
          Tup = squeeze(trapz(T2.*alfa2.*exp(-Tauupw21./UMU),1));
          Tau = mean(exp(-Taudw21./UMU));
          
          Tup = Tup';
          Tau = Tau';
%     
LST = SKINT(:,day);
load(['/Volumes/G-RAIDT/Emissivity-AMSR2/BT/2016/AMSR2_' num2str(FF1(f), '%02d') 'V_ASC_' num2str(mm(id), '%02d') '_2016.mat'])
BTEMPV = TB_V;
TBV = BTEMPV(:,day);
%
load(['/Volumes/G-RAIDT/Emissivity-AMSR2/BT/2016/AMSR2_' num2str(FF1(f), '%02d') 'H_ASC_' num2str(mm(id), '%02d') '_2016.mat'])
BTEMPH = TB_H;
TBH = BTEMPH(:,day);
%%%%%%%  Emissivity computation  %%%%%%%%%%%%%%%%%%%%%%%%
NT3 = Tup.*Tau;
DT2 = LST-Tup;
NUMV = TBV-Tup-NT3;
NUMH = TBH-Tup-NT3;
DEN = Tau.*DT2;

EMISV = NUMV./DEN;
EMISV(isnan(LST) | isnan(TBV) | isnan(WV)) = NaN;
emis_v(:,day,f) = EMISV;
%
EMISH = NUMH./DEN;
EMISH(isnan(LST) | isnan(TBH) | isnan(WV)) = NaN;
emis_h(:,day,f) = EMISH;

end;

end;
%%%%%%%%%%%%%%%%%%%% WRITING OUTPUT FILES %%%%%%%%%%%%%%%%%%%%%
save(['/Volumes/G-RAIDT/Emissivity-AMSR2/Emissivity/2016/Uncorrected_V1/AMSR2V_EMIS_ASC_' num2str(mm(id), '%02d') '_2016.mat'],'emis_v');
save(['/Volumes/G-RAIDT/Emissivity-AMSR2/Emissivity/2016/Uncorrected_V1/AMSR2H_EMIS_ASC_' num2str(mm(id), '%02d') '_2016.mat'],'emis_h');
end;
%
% %%%%%%%%%%%  Restructing for plotting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ilat cellcntr icells box flat flon dlont thismax iind jind]=textread(...
%     '/volumes/G-RAID1/Mobile2/useful/p25ancil.out','%d%d%d%d%f%f%f%f%d%d');
% cellN = load('/volumes/G-RAID1/Mobile2/Useful/LandcellN.dat');
% 
% TT=zeros(660066,1);
% TT(cellN)=emis_h(:,31,7);
% TT(TT == 0) = NaN;
% %
% mtx=zeros(1440,720);
% for i=1:1440*720
%     mtx(i)= TT(box(i));
% end;
% mtx1=mtx;
% mtx1(1:720,:)=mtx(721:1440,:);
% mtx1(721:1440,:)=mtx(1:720,:);
% %%%%%%%%%%%%%%%% PLOTTING FIGURE %%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% imagesc(flipud(mtx1'));
% colormap(jet);
% colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%