function PerformanceEmissionsPlotting(ax,a,b,y1,y2)
    
    load("E:\PAPERS\XXXX H2DDI Intake boosting\20250930 H30bTDC\Data.mat");
    
    axes(ax);


    C = orderedcolors("gem12");
    C = [0 0 0; C];

    

    hold on;
    grid on;

    xlim([(Data{2,18}-5) Data{4,18}+5]);
    xticks((Data{2,18}):10:(Data{4,18}));
    %xlabel('First H_2 fraction of total 95% energy [%]','FontSize', 12, 'FontWeight', 'bold');
    
    %Left data  
    %yyaxis left;
    set(ax, 'ColorOrder', C);
    %ax.YColor = [0	0.447 0.741];
    %plot(ax,[(Data{1,x}),(Data{4,x})],[Data{1,y1},Data{1,y1}],LineWidth=1.5,LineStyle="--",Color=C(4,:),Marker="none");
    plot(ax,Data{a:b,18},Data{a:b,y1},LineWidth=1.5,LineStyle="-",Marker='o',Color=C(y2,:));
   
    %%Right Data
    %yyaxis right;
    %set(ax, 'ColorOrder', C);
    %ax.YColor = [0.85 0.325 0.098];
    %%plot(ax,[(Data{2,x}),(Data{end,x})],[Data{1,y2},Data{1,y2}],LineWidth=1.5,LineStyle="--",Color=C(3,:),Marker="none");
    %plot(ax,Data{a:b,18},Data{a:b,y2},LineWidth=1.5,LineStyle="-",Marker="o",Color=C(3,:));

end