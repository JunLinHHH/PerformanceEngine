function PressureAndaHRR(ax,TestDate,SetNr,i) 
    
    %% Fill the processed data folder path
    Enginedatafolder = "E:\H2DDI Engine Data";
    
    %% Automaticlly find the data file Date_Set_No_RPM_Ti_Tw
    
    rpm = 1400; % Engine speed
    Ti = 30; % Inlet temperature
    Tw = 90; % Coolent temperature
    
    Processeddatafolder_name = ['_','Set',num2str(SetNr,'%02d'),'_',num2str(rpm),'_',num2str(Ti),'_',num2str(Tw)];
    Processeddata_name = ['_','Set',num2str(SetNr,'%02d'),'_',num2str(rpm),'_',num2str(Ti),'_',num2str(Tw),'_','processed.mat'];
    Processeddata_path = fullfile(Enginedatafolder,num2str(TestDate,'%02d'),Processeddatafolder_name,Processeddata_name);
    disp(Processeddata_path);
    loadedData= load(Processeddata_path,'Stats','Thermo','Traces');

    C = orderedcolors("gem12");
    C = [0 0 0; C];
    
    %% Pressure and aHRR doubleY
    axes(ax);
    
    hold on;
    grid on;
    title('Title','FontSize', 14, 'FontWeight', 'bold');

    xlim([-20 40])
    xticks(-20:10:40)
    xlabel('Crank Angle [°CA aTDC]','FontSize', 14, 'FontWeight', 'bold');

    % Pressure  
    yyaxis left;
    
    ylim([-10.5 10])
    yticks(-10:1:10)
    ylabel('Pressure [MPa]', 'FontSize', 14, 'FontWeight', 'bold','Units','data','Position',[-26,6,0]);
    set(ax,'YtickLabel',{'','','','','','','','','','','','','2','','4','','6','','8','','10'});
    ax.YColor = [0	0.447 0.741];

    Pressure = plot(ax,loadedData.Traces.CA , (loadedData.Traces.Pmean)/10,LineWidth=1.5,LineStyle="-",Color=C(i,:),Marker="none");   
    
    % aHRR
    yyaxis right;
    
    ylim([-12.5 500]);
    yticks(0:25:500)
    ylabel('aHRR [J/°CA]', 'FontSize', 14, 'FontWeight', 'bold','Units','data','Position',[48,150,0]);
    set(ax,'Ytick',0:25:500,'YtickLabel',{'0','','50','','100','','150','','200','','250','','300','','','','','','','',''});
    ax.YColor = [0.85 0.325 0.098];
    
    loadedData.Thermo.HRRfilt = bandstop(loadedData.Thermo.aHRR,[0.4 0.7],10); % Filter
    loadedData.Thermo.HRRfilt = lowpass(loadedData.Thermo.HRRfilt,[4],10); % Fliter

    aHRR = plot(ax,loadedData.Traces.CA(1:end-1) , loadedData.Thermo.HRRfilt,LineWidth=1.5,LineStyle="-",Color=C(i,:),Marker="none"); %Origional
    
    % Hide aHRR legend
    set(get(get(aHRR,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(gca,'FontSize',14,'FontWeight','bold');

    DSOI = loadedData.Traces.CA(find(loadedData.Traces.InjMean>0.01,1))+2.52; % 2.52 is the delay of the actuall diesel injection = 0.3ms * 8.4CA/ms
    DEOI = loadedData.Traces.CA(find(loadedData.Traces.InjMean>0.01,1,'last'))+2.52;
    Inj_dur = DEOI - DSOI;
    
    TextMarking={'Diesel/100kPa','100kPa','110kPa','120kPa','119kPa','50/50/D','40/60/D','30/70/D','20/80/D','10/90/D'};
    
    %plot(ax,DEOI,(82-10*i),'>','Color',C(i,:),'MarkerFaceColor',C(i,:)); %末尾箭头
    %plot(ax,[DSOI DEOI],[1 1 ]*(82-12*i),'-','Color',C(i,:)); %中间横线
    %plot(ax,[DSOI DSOI],[-20,(87-12*i)],'-','Color',C(i,:)); %起始时刻的竖线
    %text(ax,-0.2,82-10*i,TextMarking{i},'HorizontalAlignment','right','FontSize',10,'Color',C(i,:));
    
    % Diesel injection timing line
    baseY   = 120;     % 第一行的 y 位置
    gapY    = 15;     % 行间距（想更稀疏就调大）
    
    %% 然后在后面的绘图语句里统一使用
    
    yLine   = baseY - gapY*(i);      % 水平线和文字的 y
    yText   = yLine + 2;               % 文字跟水平线对齐
    yVLine  = -20;                 % 竖线起点
    yVTop   = yLine + 8;           % 竖线顶端（可视需要加一点）
    
    plot(ax, DEOI,        yLine,          '>','Color',C(i,:), 'MarkerFaceColor',C(i,:));
    
    plot(ax,[DSOI DEOI], [yLine yLine],   '-', 'Color',C(i,:));
    
    plot(ax,[DSOI DSOI], [yVLine yVTop],  '-', 'Color',C(i,:));
    
    text(ax, DSOI-0.2, yText, TextMarking{i}, 'HorizontalAlignment','right','FontSize',10,'Color',C(i,:));

end

%%lgd = legend('Diesel/DSOI = -5','0/DSOI = -7','-10/DSOI = -9','-20/DSOI = -4','-40/DSOI = -7','-60/DSOI = -7','-120/DSOI = -5','-180/DSOI = -4','Fontsize',8,'Fontweight','normal','Location','east');
%title(legend,sprintf('MeOH SOI\n[°CA bTDC]'))
%title(legend,sprintf('MeOH fraction\n [%%]'))
%% Hydrogen Legend
%title(sprintf('H_2 fraction 60%%'))
%title(sprintf('H_2SOI = 150 °CA bTDC'))

%lgd = legend('Diesel/2','0/2','60/4','120/3','-40/DSOI = -7','-60/DSOI = -7','-120/DSOI = -5','-180/DSOI =-4','Fontsize',8,'Fontweight','normal','units','pixels','Location',[180 120 90 100]);
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

%arrow([1, 42.5], [5.5, 52.5], ...
%    'Length', 40, ...
%    'TipAngle', 25, ...
%    'BaseAngle', 90, ...
%    'FaceColor', 'none', ...
%    'EdgeColor', 'k', ...
%    'LineWidth', 2,'width',20)