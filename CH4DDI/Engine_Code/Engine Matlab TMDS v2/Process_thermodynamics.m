function handles=Process_thermodynamics(handles)
Kappa = 1.40;

DV = 0.005890/6; %displacement volume [m^3]
CR = 17.4; %compression ratio
CV = DV/(CR - 1); %clearance volume [m^3]

Pmean = handles.traces.Pmean;
CA = handles.traces.CA;
Vol = handles.traces.Vol;

%% Estimate the pressure sensor offset
CA_Range = find(CA>-143,1):find(CA>-110,1);
fun  = @(P0) sum((Pmean(CA_Range)' - P0(1) - P0(2)*((CV+DV)./Vol(CA_Range)).^Kappa).^2);
P0_guess = [-1, 1];
P0 = fminsearch(fun,P0_guess);
P0(3)=Kappa;

%Alternative estimation if fit gives unreasonable values, set the pressure
%in the BDC to 1
% if P0(2)>1.15 || P0(2)<0.9
    CA_Range = find(CA>-210,1):find(CA>-180,1);
    P0(1) = mean(Pmean(CA_Range))-0.1;
    Pmean=Pmean-P0(1)+0.9;

    %%
    [~, idx_TDC] = min(abs(CA)); 

    % 计算平均 TDC 处的压力
    P_TDC_mean = Pmean(idx_TDC, :); % 取出所有循环mean的 CA=0° 处的压力
    

% else
%     Pmean = Pmean-P0(1)-0.2;
% end

handles.traces.Pmean = Pmean;
handles.traces.P_TDC_mean = P_TDC_mean; % 存储平均 TDC 处压力
Kappa = 1.35;

% Calculate the aHRR
aHRR = (Kappa/(Kappa-1)*Pmean(1:end-1)'.*diff(Vol) + 1/(Kappa-1)*Vol(1:end-1).*diff(Pmean'));
aHRR = aHRR *1e5/(CA(2)-CA(1));
cumHRR = cumsum(aHRR)*(CA(2)-CA(1));
cumHRR = cumHRR - cumHRR(find(CA>-50,1));
CA_Range = find(CA>-50,1):find(CA>-30,1);
minCum = min(cumHRR(CA_Range));
CA_Range = find(CA>0,1):find(CA>180,1);
maxCum = max(cumHRR(CA_Range));

SumHRR = maxCum-minCum;
CA_Range = find(CA>-50,1):find(CA>150,1);
CA10 =CA(find(CA>-50,1)+find(cumHRR(CA_Range)>minCum+0.1*SumHRR,1,'first'));
CA50 =CA(find(CA>-50,1)+find(cumHRR(CA_Range)>minCum+0.5*SumHRR,1,'first'));
CA90 =CA(find(CA>-50,1)+find(cumHRR(CA_Range)>minCum+0.8*SumHRR,1,'first'));
BurnDur = CA90-CA10;
IgnDelay = CA10 - CA(find(handles.traces.InjMean>0.5,1,'first'));

% Store values into handles
handles.thermo.aHRR = aHRR;
handles.thermo.cumHRR = cumHRR;
handles.thermo.totHRR = SumHRR;
handles.thermo.CA10 = CA10;
handles.thermo.CA50 = CA50;
handles.thermo.CA90 = CA90;
handles.thermo.BurnDur = BurnDur;
handles.thermo.IgnDelay = IgnDelay;
end