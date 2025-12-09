% clear all


Tw = 90; %C
Ti = 30; %C
rpm = 1400; %RPM
Pinj = 1000; %ba
Ptinj =  -0; %aTDC
Pdur = 0; %microseconds
Mtinj = -5; %aTDC
Mdur = 650; %microseconds


HC_M = 12*6+14*1;      % Hexane Molecular weight
Air_M = 28.84;         % Air Molecular weight
NO_M = 14*1+16*1;
NO2_M = 14*1+16*2;
CO_M = 12*1+16*1;
CO2_M = 12*1+16*2;
afr_vol = 981.67e-6*(rpm/2/60) ;   %m3 per sec



storename = strcat('E:\PAPERS\XXXX H2DDI Intake boosting\20251030 H2SOI28-24 - update HSOI 20\',num2str(rpm),'_',num2str(Tw),'_',num2str(Ti),'_',num2str(Pinj),'_',num2str(Ptinj),'_',num2str(Pdur),'_',num2str(Mtinj),'_',num2str(Mdur),'_Torque_variation','.mat');

% mpi=14.817; % mg/injection
CV_Diesel = 43.4; % MJ/kg;
% CV_Ethanol= 26.87 % MJ/kg;
% CV_Gasoline= 44.4 % MJ/kg
DV = 0.005890/6;  %m^3
rho_air = 1.184;       % density of air at 25 degree C
M_air = 28.97; % Molecular weight of air
M_NO = 30;    % Molecular weight of NO
M_NO2 = 46;    % Molecular weight of NO2
load("E:\PAPERS\XXXX H2DDI Intake boosting\20251030 H2SOI28-24 - update HSOI 20\Emissions&Eff.mat")
IAT = GCI_Tor(:,1);
GCI_fuel = GCI_Tor(:,2);    %Fuel mass in [mg]
Pinjection = GCI_Tor(:,3); %Injection pressure
Ptiming = GCI_Tor(:,4);    %Primary injection timing
Pduration =  GCI_Tor(:,5); %Primary injection duration
Mtiming = GCI_Tor(:,6);    %Main injection timing 
Mduration = GCI_Tor(:,7);  %Main injection duration
BrakeTor = GCI_Tor(:,8);
IMEP_Kpa = GCI_Tor(:,9);
NOx_ppm = GCI_Tor(:,10);        %ppm
NO_ppm = GCI_Tor(:,11);         %ppm
CO_vol = GCI_Tor(:,12);
HC_ppm = GCI_Tor(:,13);        %ppm      
CO2_vol = GCI_Tor(:,14);      
Opacity = GCI_Tor(:,15);
CE = GCI_Tor(:,16);
MAP = GCI_Tor(:,17);

R = 287.05;
P = MAP*1e5;                 % Intake pressure MAP
T = 25 + 273.15;
rho = P / (R * T);
afr_mass = afr_vol*rho*1000;           % rho density g/m3 at 25 celus for each MAP
afr_mol = afr_mass/28.96;  

NO2_ppm = NOx_ppm-NO_ppm;

IP = IMEP_Kpa*DV*(rpm/120); %KW % step1
BP = 2*pi()*rpm*BrakeTor/60/1000; %kW
HC_gmkwH = ((HC_ppm/1e6).*afr_mol)*HC_M./IP*3600;
CO_gmkwH = ((CO_vol/100).*afr_mol)*CO_M./IP*3600;
CO2_gmkwH = ((CO2_vol/100).*afr_mol)*CO2_M./IP*3600;

NOx_gmkWH = (((NO_ppm/1e6).*afr_mol)*(NO_M)+((NO2_ppm/1e6).*afr_mol)*(NO2_M))./IP*3600;
NO2_gmKWH = ((NO2_ppm/1e6).*afr_mol)*(NO2_M)./IP*3600;
NO_gmKWH = ((NO_ppm/1e6).*afr_mol)*(NO_M)./IP*3600;

fuel = GCI_fuel;
imep = IMEP_Kpa;

bmep = BP*1000/(DV*rpm/2/60)/1000;        % kPa  
fmep = imep-bmep;
BTheff = BP./((GCI_fuel/1e6)*rpm/120*CV_Diesel*1e3)*100;  %%
ITheff = IP./((GCI_fuel/1e6)*rpm/120*CV_Diesel*1e3)*100;  %%

BSFC = (fuel/1000*rpm/2/60)./BP*3600;  % g/kW-H
ISFC = (fuel/1000*rpm/2/60)./IP*3600;  % g/kW-H

save (storename, 'IAT','Pinjection','Ptiming','Pduration','Mtiming','Mduration','BrakeTor', 'NOx_ppm', 'NO2_gmKWH','NO_ppm','NO_gmKWH','NO2_ppm','IP', 'BP','bmep','fmep','BTheff','ITheff','ISFC','BSFC','HC_ppm','CO_vol','CO2_gmkwH','Opacity','HC_gmkwH', 'NOx_gmkWH', 'CO_gmkwH','CE');
% save (storename, 'IAT','Pinjection','Ptiming','Pduration','Mtiming','Mduration','BrakeTor','IP', 'BP','bmep','fmep','BTheff','ITheff','ISFC','BSFC');
