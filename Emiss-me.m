%%%% Computation of emissivity at AMSR-2 frequencies after computing 
%%%% atmospheric extinction, upwelling and downwelling Tbs from
%%%%   AIRS air temperaure (K) profile and TWV (kg/m^2) 
%%%%        NOTE:  DAY--> ASCENDING   &  NIGHT-->DESCENDING
%%%%              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
%
FF = [6.925; 7.30; 10.65; 18.7; 23.8; 36.5; 89.0];
FF1 = [06; 07; 10; 18; 23; 36; 89];

UMU = cosd(55);  %AMSR-2 incidence angle
%

load(['/Users/student/Documents/NASA Project/TotH2OVap_A_Jan2014.mat']);
TWV = LandTWV./10;  %%Conveting kg/m^2 to cm
%
load(['/Users/student/Documents/NASA Project/airt_Jan2014.mat']);
AIRT = LANDAirT;
%
load(['/Users/student/Documents/NASA Project/LST_Jan2014.mat']);
SKINT = LANDT;
%
for day=1:31
    
    
WV = TWV(:,day);
        
Temp_Profile = AIRT(:,:,1);
 
for f=1:7
    freq=FF(f)
   

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
load(['/Users/student/Documents/NASA Project/' num2str(FF1(f), '%02d') 'V_ASC_Jan2014.mat'])
BTEMPV = TB_V;
TBV = BTEMPV(:,day);
%
load(['/Users/student/Documents/NASA Project/' num2str(FF1(f), '%02d') 'H_ASC_Jan2014.mat'])
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
%save(['/Users/student/Documents/NASA Project/AMSR2V_EMIS_ASC_JAn2014.mat'],'emis_v');
%save(['/Users/student/Documents/NASA Project/AMSR2H_EMIS_ASC_Jan2014.mat'],'emis_h');

%
% %%%%%%%%%%%  Restructing for plotting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ilat cellcntr icells box flat flon dlont thismax iind jind]=textread('p25ancil.out','%d%d%d%d%f%f%f%f%d%d');
cellN = load('LandcellN.dat');

TT=zeros(660066,1);
TT(cellN)=emis_h(:,31,7);
TT(TT == 0) = NaN;
%
mtx=zeros(1440,720);
for i=1:1440*720
    mtx(i)= TT(box(i));
end;
mtx1=mtx;
mtx1(1:720,:)=mtx(721:1440,:);
mtx1(721:1440,:)=mtx(1:720,:);
% %%%%%%%%%%%%%%%% PLOTTING FIGURE %%%%%%%%%%%%%%%%%%%%%%%%
figure;
imagesc(flipud(mtx1'));
colormap(jet);
colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%