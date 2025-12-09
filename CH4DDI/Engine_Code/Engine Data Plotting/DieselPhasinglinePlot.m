% Diesel plot

function DieselPhasinglinePlot(ax11,ax12,ax13,ax21,ax22,ax23)

        FolderPath = 'C:\Users\82050\OneDrive - UNSW\UNSW_PhD_Chris_SK\H2DDI\1400J\Processed Figures\20241201';
        
        DieselDataPath = fullfile(FolderPath,'H0Phasing.mat');

        DieselData = load(DieselDataPath);



       

        plot(ax21,[-160,10],[DieselData.DataForPlot(5),DieselData.DataForPlot(5)],linestyle="--",LineWidth=2,Color=[96 96 96]/255, HandleVisibility="off");
        
        lgd = legend('60% H_2','80% H_2','95% H_2','Fontsize',8,'Fontweight','normal','units','pixels','Location',[80 525 90 50]);

        plot(ax22,[-160,10],[DieselData.DataForPlot(6),DieselData.DataForPlot(6)],linestyle="--",LineWidth=2,Color=[96 96 96]/255, HandleVisibility="off");
        

        plot(ax23,[-160,10],[DieselData.DataForPlot(7),DieselData.DataForPlot(7)],linestyle="--",LineWidth=2,Color=[96 96 96]/255, HandleVisibility="off");
        

        plot(ax11,[-160,10],[DieselData.DataForPlot(6) - DieselData.DataForPlot(5),DieselData.DataForPlot(6) - DieselData.DataForPlot(5)],linestyle="--",LineWidth=2,Color=[96 96 96]/255, HandleVisibility="off");
       

        plot(ax12,[-160,10],[DieselData.DataForPlot(7) - DieselData.DataForPlot(6),DieselData.DataForPlot(7) - DieselData.DataForPlot(6)],linestyle="--",LineWidth=2,Color=[96 96 96]/255, HandleVisibility="off");
        

        plot(ax13,[-160,10],[DieselData.DataForPlot(8),DieselData.DataForPlot(8)],linestyle="--",LineWidth=2,Color=[96 96 96]/255, HandleVisibility="off");
              



end