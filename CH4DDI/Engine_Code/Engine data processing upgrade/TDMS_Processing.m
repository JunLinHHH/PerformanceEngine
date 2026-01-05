% Ver 1.0
% Current MCCDAQ board Di4(Hydrogen injection signal) Di3(Diesel injector 
% signal); Di6(Encoder degree)
% Di7(Encoder original Z-pulse); Ai7(Pressure signal)
% Change the main data folder
%% Main function
function [EngineTestCase_info,TDMS_RawData,ProcessedData] = TDMS_Processing(TestDate,SetN,DSOI)

    % Engine operating condition 
    DV = 0.005890/6; %displacement volume [m^3]
    CR = 17.4; %compression ratio
    crlen = 0.2221; %connecting rod length [m] 
    crad = 0.0625; %crank radius [m] (Normally half of S)
    B = 0.100; %Bore
    TDC_shift = -67;
    CA = (0:0.1:719.9)-360;
    Lambda = crad / crlen;
    Vol = DV*(1/(CR-1) + 0.5*(1+1/Lambda*(1-sqrt(1-(Lambda*sind(CA)).^2))-cosd(CA)));
    rpm = 1400; %Engine speed
    IAT = 30;   %Intake Air temperature
    TW  = 90;   %Coolant temperature

    % Folder path - use char arrays instead of string arrays for exist() compatibility
    MainDATA_folder = char(fullfile('I:\CH4DDI\Code_development\Chirs data\20240906_H2DDI'));
    
    % Convert SetN to number if it's a string
    SetN_num = str2double(SetN);
    
    % Try both folder naming conventions (with and without leading zero)
    Targetdataset_folder_v1 = char(fullfile(MainDATA_folder,num2str(TestDate),['_Set',num2str(SetN_num),'_',num2str(rpm),'_',num2str(IAT),'_',num2str(TW)]));
    Targetdataset_folder_v2 = char(fullfile(MainDATA_folder,num2str(TestDate),['_Set',num2str(SetN_num,'%02d'),'_',num2str(rpm),'_',num2str(IAT),'_',num2str(TW)]));
    
    % Check which folder exists
    if exist(Targetdataset_folder_v1, 'dir')
        Targetdataset_folder = Targetdataset_folder_v1;
        fprintf('Found folder: %s\n', Targetdataset_folder_v1);
    elseif exist(Targetdataset_folder_v2, 'dir')
        Targetdataset_folder = Targetdataset_folder_v2;
        fprintf('Found folder: %s\n', Targetdataset_folder_v2);
    else
        error('Could not find data folder. Tried:\n  %s\n  %s', Targetdataset_folder_v1, Targetdataset_folder_v2);
    end
    
    % Create a structure for storing process
    EngineTestCase_info = struct();
    TDMS_RawData = struct();
    ProcessedData = struct();

    % Find TDMS files and save their path
    clc;
    TDMS_files = dir(fullfile(Targetdataset_folder, '*.tdms'));
    TDMS_files_numbers = length(TDMS_files);
    disp(['Selected folder: ',Targetdataset_folder]);
    disp(['Total ',num2str(TDMS_files_numbers), ' TDMS files']);
    
    % Check if any TDMS files found
    if TDMS_files_numbers == 0
        error('No TDMS files found in folder: %s', Targetdataset_folder);
    end
    
    % Save result in struct
    % Save engine operating condition to structure
    EngineTestCase_info.DV = DV;
    EngineTestCase_info.CR = CR;
    EngineTestCase_info.crlen = crlen;
    EngineTestCase_info.crad = crad;
    EngineTestCase_info.B = B;
    EngineTestCase_info.TDC_shift = TDC_shift;
    EngineTestCase_info.CA_range = CA;
    EngineTestCase_info.Lambda = Lambda;
    EngineTestCase_info.Vol_CA = Vol;
    EngineTestCase_info.Engine_Speed_rpm = rpm;
    EngineTestCase_info.T_Intake_Air = IAT;
    EngineTestCase_info.T_Coolant = TW;

    ProcessedData.CA = CA;

    % Save TDMS file paths
    TDMS_RawData.TDMSfile_paths = cell(TDMS_files_numbers,1);
    for fn = 1:TDMS_files_numbers
        TDMS_RawData.TDMSfile_paths{fn,1} = fullfile(Targetdataset_folder,TDMS_files(fn).name);
        disp(['File ',num2str(fn),': ' TDMS_files(fn).name])
    end

    % Save TDMS raw data 
    TDMS_RawData.TDMS_info = ReadTDMS(TDMS_RawData.TDMSfile_paths);

    % Save TDC positions
    ProcessedData.TDCs = Find_TDC(TDMS_RawData.TDMS_info,rpm,TDC_shift);
    [ProcessedData.mean_H_InjT,ProcessedData.mean_InjT,ProcessedData.Pmean,ProcessedData.P_single,ProcessedData.IMEP,ProcessedData.CoV,ProcessedData.P_max,ProcessedData.PRR_max] = P_single(TDMS_RawData.TDMS_info,ProcessedData.TDCs,CA,Vol,DV,DSOI);
    [ProcessedData.TDC_P_mean,ProcessedData.Pmean,ProcessedData.aHRR,ProcessedData.cumHRR,ProcessedData.SumHRR,ProcessedData.CA10,ProcessedData.CA50,ProcessedData.CA90,ProcessedData.BurnDur,ProcessedData.IgnDelay] = Process_thermodynamics(ProcessedData,DV,CR,ProcessedData.Pmean,CA,Vol);
    figure_plotting(ProcessedData.TDC_P_mean,ProcessedData.CA,Vol,ProcessedData.Pmean,ProcessedData.aHRR,ProcessedData.cumHRR,ProcessedData.mean_H_InjT,ProcessedData.mean_InjT,ProcessedData.P_max,ProcessedData.CoV,ProcessedData.CA10,ProcessedData.CA50,ProcessedData.CA90,ProcessedData.IgnDelay,ProcessedData.BurnDur,ProcessedData.IMEP,ProcessedData.PRR_max,ProcessedData.SumHRR);
    
    % Save the processed data into .mat file
    StoreProcessed_name = fullfile(Targetdataset_folder,['_',num2str(TestDate),'_','Set_',num2str(SetN,'%02d'),'_','Processed.mat']);
    save(StoreProcessed_name,'EngineTestCase_info','TDMS_RawData','ProcessedData')
end

%% Get TDMS raw data    
function TDMS_info = ReadTDMS(handles)
    TDMS_info = cell(length(handles),1);
    for i = 1 : length(handles)
        TDMS_info{i,1} = TDMS_getStruct(handles{i,1});
    end
end

%% Find TDCs
function [TDCs, Pressure_all, Injection_all] = Find_TDC(handles,rpm,TDC_shift)
    
    TDCs = cell(length(handles),1);
    Pressure_all = cell(length(handles),1);
    Injection_all = cell(length(handles),1);

    for i = 1:length(handles)
        Zpulse = double(handles{i,1}.Digital_channels.Untitled_7.data);

        % Find the TDC positions
        Zfilt = Zpulse;
        TDC_curr = find(Zfilt==1 & ([1 Zfilt(1:end-1)])==0);
        TDC_diff = mean(diff(TDC_curr));
        TDC_curr = TDC_curr - round((TDC_shift) / 720 * TDC_diff);

        %get rid of the chipped cycles
        if(TDC_curr(1)<1)
            TDC_curr = TDC_curr(2:end);
        end
        if (TDC_curr(end)>length(Zfilt))
            TDC_curr=TDC_curr(1:end-1);
        end
        TDCs{i}=TDC_curr;

        %extract the CA degree from the encoder trace
        PulsTrace = handles{i,1}.Digital_channels.Untitled_6.data;
        CA_Enc_curr{i} = cumsum(double(PulsTrace>0 & ([1 PulsTrace(1:end-1)])==0))*0.2;
        a=1; 
        b=round((60/rpm/1800)/handles{i,1}.Analog_channels.Dev0_Ai5.Props.wf_increment*5);
        b=ones(b,1)/b;
        CA_Enc_curr{i} =filtfilt(b,a,CA_Enc_curr{i});
        CA_Enc_curr{i} = mod((CA_Enc_curr{i}-mean(CA_Enc_curr{i}(TDC_curr) - (1:length(TDC_curr))*720)),720)-360;
        TDCs{i} = find(CA_Enc_curr{i}>[CA_Enc_curr{i}(2:end) CA_Enc_curr{i}(end)]);
    end
end

%% Pressure traces interception (P_single,P_mean,IMEP,COV)
function [mean_H_injection,mean_injection,mean_pressure,P_single,IMEP,CoV,P_max,PRR_max] = P_single(handles,TDC,CA,Vol,DV,DSOI)

    num_cycle = 0;
    diffTDC = [];
    for fn = 1:length(TDC)
        num_cycle = num_cycle + length(TDC{fn})-2;
        diffTDC = [diffTDC diff(TDC{fn})];
    end
    RPM = 2*60/(mean(diffTDC)*handles{fn,1}.Analog_channels.Dev0_Ai5.Props.wf_increment);
    engine_speed_rpm = RPM;
    engine_speed_rps = engine_speed_rpm/60;
    engine_speed_CAps = engine_speed_rps* 360;
    sec_per_CA = 1/engine_speed_CAps;
    sampling_rate = 200 * 1000;
    data_per_CA = sec_per_CA * sampling_rate;
    CA_per_data = 1/data_per_CA;
   
    if DSOI > 0
       inj_CA = -DSOI;
    else
       inj_CA = DSOI;
    end
   
    side_left = ceil(abs((-360 - inj_CA)/CA_per_data));
    side_right = ceil(abs((360 - inj_CA)/CA_per_data));

    %CA(end) = [];
    range_data_point = -side_left : side_right;
    range_CA = (range_data_point(1)*CA_per_data: CA_per_data :range_data_point(end) * CA_per_data) + inj_CA;
    store_pressure_by_file_number = zeros(length(CA), length(handles));
    store_injection_by_file_number = zeros(length(CA), length(handles));
    P_single = [];

    % Filter the pressure traces (7kHz cutoff)
    % f_cutoff=10000;     % Cut-off frequncy for filtering(Original)
    f_cutoff=10000;       % Test one 20240216
    fc_norm=f_cutoff/(7200*(RPM/60/2));        % Normalised cut-off frequncy for filtering
    [b,a]=butter(2,fc_norm);        % Butterworth filter
   
    for fn=1:length(handles)
        Pressure = double(handles{fn,1}.Analog_channels.Dev0_Ai7.data);
        Injection = double(handles{fn,1}.Digital_channels.Untitled_5.data);
        H_InjT = double(handles{fn,1}.Digital_channels.Untitled_4.data);
    
        diff_inj = [0,diff(Injection)];
        [~, index_inj]= find(diff_inj > 0.5);
        
        idx_pressure_segment = 1;
        
        p_voltage_segement_per_file = zeros(length(CA), length(index_inj) - 2 - 1);
        inj_segement_per_file =       zeros(length(CA), length(index_inj) - 2 - 1);    
        for ii = 2: (length(index_inj) - 1 - 1)
            index_injection = index_inj(ii);
            this_cycle_pressure = Pressure(index_injection-side_left : index_injection+side_right);
            this_cycle_pressure_filter = filtfilt(b,a,this_cycle_pressure);
            this_cycle_pressure_filter_rescale =interp1(range_CA,this_cycle_pressure_filter,CA,'linear','extrap');
            p_voltage_segement_per_file(:,idx_pressure_segment) = this_cycle_pressure_filter_rescale;
            idx_pressure_segment = idx_pressure_segment + 1;

            this_cycle_injection = Injection(index_injection-side_left : index_injection+side_right);
            this_cycle_H_injection = H_InjT(index_injection-side_left : index_injection+side_right);

            this_cycle_injection_rescale =interp1(range_CA,this_cycle_injection,CA,'linear','extrap');
            this_cycle_H_injection_rescale =interp1(range_CA,this_cycle_H_injection,CA,'linear','extrap');

            inj_segement_per_file(:,idx_pressure_segment) = this_cycle_injection_rescale;
            H_inj_segement_per_file(:,idx_pressure_segment) = this_cycle_H_injection_rescale;

        end
        pressure_segement_per_file = p_voltage_segement_per_file * 20;
        mean_pressure_by_one_file = mean(pressure_segement_per_file,2);

        mean_injection_by_one_file = mean(inj_segement_per_file,2);
        mean_H_injection_by_one_file = mean(H_inj_segement_per_file,2);

        store_pressure_by_file_number(:,fn) = mean_pressure_by_one_file;
        store_injection_by_file_number(:,fn) = mean_injection_by_one_file;
        store_H_injection_by_file_number(:,fn) = mean_H_injection_by_one_file;

        fprintf("file number %02d done\n", fn);
        P_single = [P_single, pressure_segement_per_file];
    end
    mean_pressure = mean(store_pressure_by_file_number,2);
    mean_injection = mean(store_injection_by_file_number,2);
    mean_H_injection = mean(store_H_injection_by_file_number,2);

    Pmean = mean_pressure;
    %Pmean = filtfilt(b,a,Pmean);
    IMEP = nansum(Pmean(1:end-1)' .* (diff(Vol)))/DV;
    IMEP_single = (nansum(P_single(3:end-4,:).*repmat(diff(Vol(3:end-3))',1,size(P_single,2)),1))'/DV;
    CoV = nanstd(IMEP_single)/IMEP;
    P_max = max(Pmean)-min(Pmean)+0.8;
    PRR_max = max(diff(Pmean))*10;
end

%% Thermpdynamics caulculaion
function [TDC_P_mean,Pmean,aHRR,cumHRR,SumHRR,CA10,CA50,CA90,BurnDur,IgnDelay] = Process_thermodynamics(handles,DV,CR,Pmean,CA,Vol)
    
    CV = DV/(CR - 1); %clearance volume [m^3]
    Kappa = 1.40; % specific heat capacities 

    CA_Range = find(CA>-143,1):find(CA>-110,1);
    fun  = @(P0) sum((Pmean(CA_Range)' - P0(1) - P0(2)*((CV+DV)./Vol(CA_Range)).^Kappa).^2);
    P0_guess = [-1, 1];
    P0 = fminsearch(fun,P0_guess);
    P0(3)=Kappa;
    CA_Range = find(CA>-210,1):find(CA>-180,1);
    P0(1) = mean(Pmean(CA_Range))-0.1;
    Pmean=Pmean-P0(1)+0.9;
    [~, idx_TDC] = min(abs(CA));
    TDC_P_mean = Pmean(idx_TDC, :);

    Kappa = 1.35;

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
    IgnDelay = CA10 - CA(find(handles.mean_InjT>0.5,1,'first'));

end

%% Thermpdynamics caulculaion
function figure_plotting(TDC_P_mean,CA,Vol,Pmean,aHRR,cumaHRR,HSOI,DSOI,PeakP,CoV,CA10,CA50,CA90,IgnDelay,BunDur,IMEP,PRR,TotalHRR)
    % Create a figure for plooting data
    figName = 'Engine data plotter';
    figHandle = findobj('Type', 'Figure', 'Name', figName);
    if isempty(figHandle)
        figHandle = figure('Name', figName, 'Units', 'Pixels', 'Position', [60 60 950 790]);
    else
        figure(figHandle); 
    end

    % Create multiple axes in the same figure
    % Create axes if not already existing
    axesConfig = {
        'PressureTrace', [65 450 400 300], 'Pressure trace', 'Crank Angle [°CA aTDC]', 'Pressure [bar]';
        'PVTrace', [525 450 400 300], 'P-V diagram', 'Volume[dm^3]', 'Pressure [bar]';
        'aHRRTrace', [65 65 400 300], 'aHRR trace', 'Crank Angle [°CA aTDC]', 'aHRR [J/°CA]';
    };
    % Loop through axes configurations to create or find axes
    axesHandles = struct();
    for i = 1:size(axesConfig, 1)
        tag = axesConfig{i, 1};
        position = axesConfig{i, 2};
        titleText = axesConfig{i, 3};
        xlabelText = axesConfig{i, 4};
        ylabelText = axesConfig{i, 5};

        % Find or create the axes
        ax = findobj(figHandle, 'Type', 'Axes', 'Tag', tag);
        if isempty(ax)
            ax = axes('Parent', figHandle, 'unit', 'pixels', 'position', position, ...
                'box', 'on', 'linewidth', 1.2, 'fontsize', 12, 'fontweight', 'bold', 'Tag', tag);
            title(ax, titleText);
            xlabel(ax, xlabelText);
            ylabel(ax, ylabelText);
            grid(ax, 'on');
        end
        axesHandles.(tag) = ax; % Store the handle
    end
    
    % Plotting
    set(groot, 'defaultAxesColorOrder', 'remove');
    hold(axesHandles.PressureTrace, 'on');
    yyaxis(axesHandles.PressureTrace, 'left');
    axesHandles.PressureTrace.YColor = [0 0 0];          
    axesHandles.PressureTrace.ColorOrder = [
    0	0	0
    0 0.447 0.741
    0.85 0.325 0.098
    0.929 0.694 0.125
    0.494 0.184	0.556
    0.466 0.674	0.188
    0.3010 0.745 0.933
    0.6350 0.078 0.184];
    ylim(axesHandles.PressureTrace,[0,100]);
    Pmean = sgolayfilt(Pmean,7,91);
    plot(axesHandles.PressureTrace,CA,Pmean,'LineStyle','-','Marker','none');
    yyaxis(axesHandles.PressureTrace, 'right');
    axesHandles.PressureTrace.YColor = [0 0 0];
    ylim(axesHandles.PressureTrace,[0,2])
    set(axesHandles.PressureTrace,'YtickLabel',{});
    plot(axesHandles.PressureTrace,CA,DSOI,'Color','r','LineStyle','-','Marker','none');
    plot(axesHandles.PressureTrace,CA,HSOI,'Color',[0.00,0.45,0.74],'LineStyle','-','Marker','none');
    xlim(axesHandles.PressureTrace,[-20,40]);

    hold(axesHandles.aHRRTrace, 'on');
    yyaxis(axesHandles.aHRRTrace, 'left');
    axesHandles.aHRRTrace.ColorOrder = [
    0	0	0
    0 0.447 0.741
    0.85 0.325 0.098
    0.929 0.694 0.125
    0.494 0.184	0.556
    0.466 0.674	0.188
    0.3010 0.745 0.933
    0.6350 0.078 0.184];
    ylim(axesHandles.aHRRTrace,[-50,300]);
    yticks(-50:50:300)
    aHRR = sgolayfilt(aHRR,7,91); % filter
    plot(axesHandles.aHRRTrace,CA(1:end-1),aHRR,'LineStyle','-','Marker','none');
    yyaxis(axesHandles.aHRRTrace, 'right');
    axesHandles.aHRRTrace.YColor = [0 0 0];
    ylim(axesHandles.aHRRTrace,[-200,1200]);
    yticks(-200:200:1200)
    plot(axesHandles.aHRRTrace,CA(1:end-1),cumaHRR,'Color','r','LineStyle','-','Marker','none');
    xlim(axesHandles.aHRRTrace,[-20,40]);

    hold(axesHandles.PVTrace, 'on');
    axesHandles.PVTrace.ColorOrder = [
    0	0	0
    0 0.447 0.741
    0.85 0.325 0.098
    0.929 0.694 0.125
    0.494 0.184	0.556
    0.466 0.674	0.188
    0.3010 0.745 0.933
    0.6350 0.078 0.184];
    plot(axesHandles.PVTrace,Vol*1000,Pmean,'LineStyle','-','Marker','none');

    % Insert combustion info
    persistent combustionInfoHandle;

    Combustion_info = sprintf([ ...
    'TDC P [bar]: %.2f\n',...
    'Peak P [bar]: %.1f\n', ...
    'CoV [%%]: %.2f\n', ...
    'CA10 [°]: %.1f\n', ...
    'CA50 [°]: %.1f\n', ...
    'CA90 [°]: %.1f\n', ...
    'Ign Delay [°]: %.1f\n', ...
    'Burn Dur [°]: %.1f\n', ...
    'IMEP [bar]: %.3f\n', ...
    'P Rise Rate [bar]: %.1f\n', ...
    'Total HRR [J]: %.0f\n'], ...
    TDC_P_mean,PeakP, CoV * 100, CA10, CA50, CA90, IgnDelay, BunDur, IMEP, PRR, TotalHRR);

    % Delete old annotation if it exists
    if ~isempty(combustionInfoHandle) && isvalid(combustionInfoHandle)
        delete(combustionInfoHandle); % 删除旧的 annotation
    end

    combustionInfoHandle = annotation(figHandle, 'textbox', [0.6, 0.13, 0.2, 0.3], ...
    'String', Combustion_info, 'FontSize', 12, 'FontWeight', 'bold', ...
    'BackgroundColor', 'white', 'EdgeColor', 'black', 'HorizontalAlignment', 'right');

end