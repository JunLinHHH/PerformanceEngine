% Processed TDMS data ploter
% Presure&aHRR in one figure;
% combustion phasing: CA10 CA50 and CA90 with each of burn duration of
% CA10-50, CA50-90 and CA10-90 in one figure;
% Performance: TDC pressure, Peak pressure, IMEP, Efficiency, Ignition
% delay, CoV of IMEP in one figure;
% Emissions: CO2 and NOx in one figure;

function TDMS_plotter(ax,TestDate,SetNr) 
    
    %% Fill the processed data folder path
    Enginedatafolder = "C:\Users\82050\OneDrive - UNSW\Desktop\H2DDI";
    
    %% Automaticlly find the data file Date_Set_No_RPM_Ti_Tw
    
    rpm = 1400; % Engine speed
    Ti = 30; % Inlet temperature
    Tw = 90; % Coolent temperature
    
    Processeddatafolder_name = ['_','Set',num2str(SetNr,'%02d'),'_',num2str(rpm),'_',num2str(Ti),'_',num2str(Tw)];
    Processeddata_name = ['_',num2str(TestDate),'_','Set','_',num2str(SetNr,'%02d'),'_','Processed.mat'];
    Processeddata_path = fullfile(Enginedatafolder,num2str(TestDate),Processeddatafolder_name,Processeddata_name);
    disp(Processeddata_path);
    loadedData= load(Processeddata_path,'EngineTestCase_info','ProcessedData','TDMS_RawData');

    %% Pressure and aHRR doubleY
    axes(ax);
    hold on;
    grid on;
    ax_colororder = [
        0	0	0
        0 0.447 0.741
        0.85 0.325 0.098
        0.929 0.694 0.125
        0.494 0.184	0.556
        0.466 0.674	0.188
        0.3010 0.745 0.933
        0.6350 0.078 0.184];
    
    title('Title','FontSize', 12, 'FontWeight', 'bold')
    xlim([-10 40])
    xticks(-10:10:40)
    xlabel('Crank Angle [°CA aTDC]','FontSize', 12, 'FontWeight', 'bold');

    % Pressure  
    yyaxis left;
    set(ax,'ColorOrder', ax_colororder);
    ylim([-10.5 10])
    yticks(-10:1:10)
    ylabel('Pressure [MPa]', 'FontSize', 12, 'FontWeight', 'bold','Units','data','Position',[-15,6,0]);
    set(ax,'YtickLabel',{'','','','','','','','','','','','','2','','4','','6','','8','','10'});
    ax.YColor = [0	0.447 0.741];
    loadedData.ProcessedData.Pmean = sgolayfilt(loadedData.ProcessedData.Pmean,7,91);
    Pressure = plot(ax,loadedData.ProcessedData.CA , (loadedData.ProcessedData.Pmean)/10,LineWidth=1.5);   
    
    % aHRR
    yyaxis right;
    set(ax, 'ColorOrder', ax_colororder);
    ylim([-12.5 500]);
    yticks(0:25:300)
    ylabel('aHRR [J/°CA]', 'FontSize', 12, 'FontWeight', 'bold','Units','data','Position',[46,150,0]);
    set(ax,'Ytick',0:25:300,'YtickLabel',{'0','','50','','100','','150','','200','','250','','300','','350','',''});
    ax.YColor = [0.85 0.325 0.098];
    loadedData.ProcessedData.aHRR = sgolayfilt(loadedData.ProcessedData.aHRR,7,91);
    aHRR = plot(ax,loadedData.ProcessedData.CA(1:end-1) , loadedData.ProcessedData.aHRR,LineWidth=1.5); %Origional
    
    % Hide aHRR legend
    set(get(get(aHRR,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    

end

%%lgd = legend('Diesel/DSOI = -5','0/DSOI = -7','-10/DSOI = -9','-20/DSOI = -4','-40/DSOI = -7','-60/DSOI = -7','-120/DSOI = -5','-180/DSOI = -4','Fontsize',8,'Fontweight','normal','Location','east');
%title(legend,sprintf('MeOH SOI\n[°CA bTDC]'))
%title(legend,sprintf('MeOH fraction\n [%%]'))
%% Hydrogen Legend
%title(sprintf('H_2 fraction 60%%'))
%title(sprintf('H_2SOI = 150 °CA bTDC'))

%lgd = legend('Diesel/2','0/2','60/4','120/3','-40/DSOI = -7','-60/DSOI = -7','-120/DSOI = -5','-180/DSOI =-4','Fontsize',8,'Fontweight','normal','units','pixels','Location',[230 120 110 155]);
%lgd = legend('Diesel/2','60%/2','80%/4','95%/3','-40/DSOI = -7','-60/DSOI = -7','-120/DSOI = -5','-180/DSOI =-4','Fontsize',8,'Fontweight','normal','units','pixels','Location',[180 120 90 70]);

%title(legend,sprintf('MH_2SOI/DSOI\n[°CA bTDC]'))
%title(legend,sprintf('H_2fraction\n DSOI [°CA bTDC]'))
% legend('unit','pixels','Location', [180 65 100 100]);
%%
%set(gca, 'Position', originalPosition);
%set(gca, 'XTickLabel', {'', '', '','','',''});
%set(gca, 'XTickLabel', {'', '', '','','',''});
%set(gca, 'XTickLabel', {'', '', '','','',''});
%yyaxis right
%set(gca, 'YTickLabel', {'', '', '','','',''});

