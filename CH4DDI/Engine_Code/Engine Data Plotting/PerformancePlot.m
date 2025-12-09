%PhasingPlot
function PerformancePlot(ax11,ax12,ax13,ax21,ax22,ax23,Fraction)
        
        FolderPath = 'C:\Users\82050\OneDrive - UNSW\Desktop\H2DDI\20241209\Figures\Fixed first injection150';
        FileName = ['H',num2str(Fraction),'Phasing.mat'];
        FilePath = fullfile(FolderPath,FileName);
        
        loadedans = load(FilePath);
        
        disp(loadedans.DataForPlot);
        %% IMEP
        axes(ax11);
        grid on;
        hold on;
        yyaxis left
        plot(-loadedans.DataForPlot(:,2),loadedans.DataForPlot(:,3)*100,"Marker","o",MarkerSize=10,LineWidth=1.5);
        

        xlim([-160,10]);
        xticks(-180:30:0)
        xlabel('H_2SOI [°CA aTDC]','FontSize',12,'FontWeight','bold')

        ylim([300 , 900]);
        yticks(300:100:900);
        ylabel(sprintf('Net IMEP [kPa]'),'FontSize',12,'FontWeight','bold');
        
        %legend('10% MeOH','30% MeOH','50% MeOH','70% MeOH','90% MeOH','Diesel','Fontsize',10,'Fontweight','normal','Location','NorthWest');
        
        
        
        %% COV
        axes(ax13);
        grid on;
        hold on;

        plot(-loadedans.DataForPlot(:,2),loadedans.DataForPlot(:,4)*100,"Marker","o",MarkerSize=10,LineWidth=1.5);

        xlim([-160,10]);
        xticks(-180:30:0)
        xlabel('H_2SOI [°CA aTDC]','FontSize',12,'FontWeight','bold')

        ylim([0 , 6]);
        yticks(0:1:6);
        ylabel(sprintf('CoV of IMEP [%%]'),'FontSize',12,'FontWeight','bold');
        
        %legend('10% MeOH','30% MeOH','50% MeOH','70% MeOH','90% MeOH','Fontsize',10,'Fontweight','normal','Location','NorthWest');
        

       
        %% IG
        axes(ax22);
        grid on;
        hold on;

        plot(-loadedans.DataForPlot(:,2),loadedans.DataForPlot(:,9) ,"Marker","o",MarkerSize=10,LineWidth=1.5);

        xlim([-160,10]);
        xticks(-180:30:0)
        xlabel('H_2SOI [°CA aTDC]','FontSize',12,'FontWeight','bold')

        ylim([0 , 30]);
        yticks(0:5:30);
        ylabel('Ignition Delay [°CA]','FontSize',12,'FontWeight','bold');
        
        %legend('10% MeOH','30% MeOH','50% MeOH','70% MeOH','90% MeOH','Fontsize',10,'Fontweight','normal','Location','NorthWest');
        %legend('50% MeOH','70% MeOH','90% MeOH','Fontsize',10,'Fontweight','normal','Location','NorthWest');

        %%
        axes(ax23);
        grid on;
        hold on;        
        
        plot(-loadedans.DataForPlot(:,2),loadedans.DataForPlot(:,11)/10 ,"Marker","o",MarkerSize=10,LineWidth=1.5);

        xlim([-160,10]);
        xticks(-180:30:0)
        xlabel('H_2SOI [°CA aTDC]','FontSize',12,'FontWeight','bold')

        ylim([4 , 10]);
        yticks(4:1:10);
        ylabel('Peak Pressure [MPa]','FontSize',12,'FontWeight','bold');
        
        %legend('10% MeOH','30% MeOH','50% MeOH','70% MeOH','90% MeOH','Fontsize',10,'Fontweight','normal','Location','NorthWest');
        %legend('50% MeOH','70% MeOH','90% MeOH','Fontsize',10,'Fontweight','normal','Location','NorthWest');
        
        

        %% TDCP
        axes(ax21);
        grid on;
        hold on;
        
        plot(-loadedans.DataForPlot(:,2),loadedans.DataForPlot(:,12) ,"Marker","o",MarkerSize=10,LineWidth=1.5);

        xlim([-160,10]);
        xticks(-180:30:0)
        xlabel('H_2SOI [°CA aTDC]','FontSize',12,'FontWeight','bold')

        ylim([40 , 70]);
        yticks(40:5:70);
        ylabel('TDC Pressure [bar]','FontSize',12,'FontWeight','bold');
        
        %legend('10% MeOH','30% MeOH','50% MeOH','70% MeOH','90% MeOH','Fontsize',10,'Fontweight','normal','Location','NorthWest');
        %legend('50% MeOH','70% MeOH','90% MeOH','Fontsize',10,'Fontweight','normal','Location','NorthWest');

end

