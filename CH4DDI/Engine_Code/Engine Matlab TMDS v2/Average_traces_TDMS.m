function handles = Average_traces_TDMS(handles)

%% Engine Data
DV = 0.005890/6; %displacement volume [m^3]
CR = 17.4; %compression ratio
crlen = 0.2221; %connecting rod length [m] 
crad = 0.0625; %crank radius [m]
B = 0.100;
TDC_shift =str2double(handles.TDCshift.String);

CA = (0:0.1:719.9)-360;
Lambda = crad / crlen;
Vol = DV*(1/(CR-1) + 0.5*(1+1/Lambda*(1-sqrt(1-(Lambda*sind(CA)).^2))-cosd(CA)));

%% Separate cycles and average
for fn=1:length(handles.raw_data)
    Zpulse = double(handles.raw_data(fn).Digital_channels.Untitled_7.data);
    % Find the TDC positions
    Zfilt = Zpulse;
    TDCs = find(Zfilt==1 & ([1 Zfilt(1:end-1)])==0);
    TDC_diff = mean(diff(TDCs));
    TDCs = TDCs - round((TDC_shift)/720*TDC_diff);

    %get rid of the chipped cycles
    if(TDCs(1)<1)
        TDCs = TDCs(2:end);
    end
    if (TDCs(end)>length(Zfilt))
        TDCs=TDCs(1:end-1);
    end
    TDC{fn}=TDCs;
    
    %extract the CA degree from the encoder trace
    PulsTrace = handles.raw_data(fn).Digital_channels.Untitled_6.data;
    CA_Enc{fn} = cumsum(double(PulsTrace>0 & ([1 PulsTrace(1:end-1)])==0))*0.2;
    a=1; b=round((60/2000/1800)/handles.raw_data(fn).Analog_channels.Dev0_Ai5.Props.wf_increment*5);
    b=ones(b,1)/b;
    CA_Enc{fn} =filtfilt(b,a,CA_Enc{fn});
    CA_Enc{fn} = mod((CA_Enc{fn}-mean(CA_Enc{fn}(TDCs) - (1:length(TDCs))*720)),720)-360;
    TDC{fn} = find(CA_Enc{fn}>[CA_Enc{fn}(2:end) CA_Enc{fn}(end)]);
end

%Initialize the data-storage
Pmean = zeros(size(CA));
InjMean = zeros(size(CA));
num_cycle = 0;
diffTDC = [];

for fn=1:length(TDC)
    num_cycle = num_cycle + length(TDC{fn})-2;
    diffTDC = [diffTDC diff(TDC{fn})];
end
RPM = 2*60/(mean(diffTDC)*handles.raw_data(fn).Analog_channels.Dev0_Ai5.Props.wf_increment);
P_single = zeros(length(CA),num_cycle);
Inj_single = zeros(length(CA),num_cycle);

%average the traces
cycle_count=0;
for fn=1:length(TDC)
    Pressure = double(handles.raw_data(fn).Analog_channels.Dev0_Ai5.data);
    Injection = double(handles.raw_data(fn).Digital_channels.Untitled_5.data);    
    for i=1:length(TDC{fn})-2
        cycle_count = cycle_count +1;
        cycle_points = (TDC{fn}(i)+1):1:(mean(TDC{fn}(i+1)));
        ptr_single = interp1(CA_Enc{fn}(cycle_points),Pressure(cycle_points),CA,'linear','extrap');
        inj_single =interp1(CA_Enc{fn}(cycle_points),Injection(cycle_points),CA,'linear','extrap');
        P_single(:,cycle_count)=ptr_single;
        Inj_single(:,cycle_count)=inj_single;
   
    end
end
InjMean = nanmean(Inj_single,2);
SOImean = find(Inj_single>0.5,1,'first');
% for i=1:size(P_single,2)
%     SOI = find(Inj_single(:,i)>0.5,1,'first');
%     Inj_single(:,i) = circshift(Inj_single(:,i),-(SOI-SOImean));
%     P_single(:,i) = circshift(P_single(:,i),-(SOI-SOImean));
% end
Pmean = nanmean(P_single,2)*str2double(handles.Psens.String);
InjMean = nanmean(Inj_single,2);
P_single = P_single*str2double(handles.Psens.String);


%% Filter the pressure traces (7kHz cutoff)
%f_cutoff=10000;     % Cut-off frequncy for filtering(Original)
f_cutoff=7000; 
fc_norm=f_cutoff/(7200*(RPM/60/2));        % Normalised cut-off frequncy for filtering
[b,a]=butter(2,fc_norm);        % Butterworth filter

Pmean = filtfilt(b,a,Pmean);
P_single = filtfilt(b,a,P_single);

%% Calculate statistics
IMEP = nansum(Pmean(1:end-1)' .* (diff(Vol)))/DV;

IMEP_single = (nansum(P_single(3:end-4,:).*repmat(diff(Vol(3:end-3))',1,size(P_single,2)),1))'/DV;
CoV = nanstd(IMEP_single)/IMEP;

P_max = max(Pmean)-min(Pmean)+0.8;
PRR_max = max(diff(Pmean))*10;
[~, idx_TDC] = min(abs(CA));     
TDC_P_mean = Pmean(idx_TDC,:) -min(Pmean) + 0.8 ;
%% Store the data into handles
handles.traces.Pmean = Pmean;
handles.traces.TDC_P_mean = TDC_P_mean;
handles.traces.P_single = P_single;
handles.traces.Vol = Vol;
handles.traces.CA = CA;
handles.traces.InjMean = InjMean;




handles.stats.RPM = RPM;
handles.stats.IMEP = IMEP;
handles.stats.CoV = CoV;
handles.stats.P_max = P_max;
handles.stats.PRR_max = PRR_max;

end