clear all;

%% Load the data (the path has to be set where the data is stored)
%data{1} as Diesel reference case
data{1}= load(fullfile("D:\H2DDI Data\20250618\_Set03_1400_30_90\_Set03_1400_30_90_processed.mat"),'Thermo','Traces');
data{2}= load(fullfile("D:\H2DDI Data\20250618\_Set13_1400_30_90\_Set13_1400_30_90_processed.mat"),'Thermo','Traces');
data{3}= load(fullfile("D:\H2DDI Data\20250618\_Set15_1400_30_90\_Set15_1400_30_90_processed.mat"),'Thermo','Traces');
data{4}= load(fullfile("D:\H2DDI Data\20250618\_Set17_1400_30_90\_Set17_1400_30_90_processed.mat"),'Thermo','Traces');
data{5}= load(fullfile("D:\H2DDI Data\20250618\_Set19_1400_30_90\_Set19_1400_30_90_processed.mat"),'Thermo','Traces');
data{6}= load(fullfile("D:\H2DDI Data\20250618\_Set21_1400_30_90\_Set21_1400_30_90_processed.mat"),'Thermo','Traces');
%data{7}= load(fullfile("D:\H2DDI Data\20250618\_Set23_1400_30_90\_Set23_1400_30_90_processed.mat"),'Thermo','Traces');
%data{8}= load(fullfile("D:\H2DDI Data\20250618\_Set11_1400_30_90\_Set11_1400_30_90_processed.mat"),'Thermo','Traces');
%data{9}= load(fullfile("D:\H2DDI Data\20250618\_Set06_1400_30_90\_Set06_1400_30_90_processed.mat"),'Thermo','Traces');
%data{10}= load(fullfile("D:\H2DDI Data\20250618\_Set08_1400_30_90\_Set08_1400_30_90_processed.mat"),'Thermo','Traces');

Colour1 =  [0.50 0.50 0.50]; %Diesel
Colour2 =  [    0         0.4470    0.7410]; %10/90
Colour3 =  [    0.8500    0.3250    0.0980]; %20/80
Colour4 =  [    0.9290    0.6940    0.1250]; %30/70
Colour5 =  [    0.4940    0.1840    0.5560]; %40/60
Colour6 =  [    0.4660    0.6740    0.1880]; %50/50
Colour7 =  [    0.3010    0.7450    0.9330]; %60/40
Colour8 =  [    0.6350    0.0780    0.1840]; %70/30
Colour9 =  [    1.0000    0.8390    0.0390]; %80/20
Colour10 = [    0.3960    0.5090    0.9920]; %90/10


%loading the data and detecting the SOI
for i=1:length(data)
    P_off = mean(data{i}.Traces.Pmean(1800:2000));
    data{i}.Traces.Pmean = data{i}.Traces.Pmean + 1 - P_off;
    SOI(i) = data{i}.Traces.CA(find(data{i}.Traces.InjMean>0.5,1))+2.52;%4CA deg is added to account for the delay in injector
    
    %filter the HRR
    data{i}.Thermo.HRRfilt = bandstop(data{i}.Thermo.aHRR,[0.4 0.7],10);
    data{i}.Thermo.HRRfilt = lowpass(data{i}.Thermo.HRRfilt,[4],10);
end

%Set the duration of diesel injection (in CA)
DOI = [5.46, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68];
%Calculate the end of H2 injection (from SOI of H2 and H2 injection
%duration)
% EOIh2 = -10 + [0 500 1000 1500 2000 2500]/83;

%% Plot the double YY plot
% setting the axes to have universal plot regardless of the platform
Fsize = 12;
figure('Unit','Pixel','Position',[100 100 600 900]);
ax1 = axes('units','pixels','Position',[60 70 280 400]);hold on; grid on; box on;

%plotting the left axis with pressure
yyaxis left
p(1) = plot(data{1}.Traces.CA,data{1}.Traces.Pmean/10,'-','COlor', Colour1,'LIneWidth',1.5);
p(2) = plot(data{2}.Traces.CA,data{2}.Traces.Pmean/10,'-','COlor', Colour2, 'LIneWidth',1.5);
p(3) = plot(data{3}.Traces.CA,data{3}.Traces.Pmean/10,'-','COlor', Colour3, 'LIneWidth',1.5);
p(4) = plot(data{4}.Traces.CA,data{4}.Traces.Pmean/10,'-','COlor', Colour4,'LIneWidth',1.5);
p(5) = plot(data{5}.Traces.CA,data{5}.Traces.Pmean/10,'-','COlor', Colour5,'LIneWidth',1.5);
p(6) = plot(data{6}.Traces.CA,data{6}.Traces.Pmean/10,'-','COlor', Colour6,'LIneWidth',1.5);
%p(7) = plot(data{7}.Traces.CA,data{7}.Traces.Pmean/10,'-','COlor', Colour7,'LIneWidth',1.5);
%p(8) = plot(data{8}.Traces.CA,data{8}.Traces.Pmean/10,'-','COlor', Colour8,'LIneWidth',1.5);
%p(9) = plot(data{9}.Traces.CA,data{9}.Traces.Pmean/10,'-','COlor', Colour9,'LIneWidth',1.5);
%p(10) = plot(data{10}.Traces.CA,data{10}.Traces.Pmean/10,'-','COlor', Colour10,'LIneWidth',1.5);

%setting the Y limits and ticks to match with the HRR plot
xlim([-10 40]);ylim([-10.5 10]); set(gca,'Xtick',-10:10:40,'Ytick',-10:1:11,'YtickLabel',{'','','','','','','','','','','','','2','','4','','6','','8','','10'});
set(gca,'FOntSize',Fsize,'FontWeight','bold');
set(get(gca,'Xlabel'),'String','Crank Angle [°CA aTDC]');
set(get(gca,'Ylabel'),'String','Pressure [MPa]');
set(get(gca,'Ylabel'),'Position',get(get(gca,'Ylabel'),'Position')+[0 5 0]);

% %Plot the Hydrogen injection timing indicators
% TextMarking={'','20% H_2','30% H_2','40% H_2','50% H_2'};
% for i=2:5
%     plot(EOIh2(i)-0.6,(86-4*(5-i)),'>','Color',get(p(i),'Color'),'MarkerFaceColor',get(p(i),'Color'));
%     plot([-50 EOIh2(i)],[1 1 ]*(86-4*(5-i)),'-','Color',get(p(i),'Color'));
%     plot([EOIh2(i) EOIh2(i)],[90-4*(5-i), 82-4*(5-i)],'-','Color',get(p(i),'Color'));
%     text(EOIh2(i)+0.5,86-4*(5-i),TextMarking{i},'HorizontalAlignment','left','FontSize',Fsize-6,'Color',get(p(i),'Color'));
% end

%plot the AHRR and other details
yyaxis right
pl(1) = plot(data{1}.Traces.CA(1:end-1),data{1}.Thermo.HRRfilt,'-','COlor', Colour1, 'LIneWidth',1.5);
pl(2) = plot(data{2}.Traces.CA(1:end-1),data{2}.Thermo.HRRfilt,'-','COlor', Colour2, 'LIneWidth',1.5);
pl(3) = plot(data{3}.Traces.CA(1:end-1),data{3}.Thermo.HRRfilt,'-','COlor', Colour3, 'LIneWidth',1.5);
pl(4) = plot(data{4}.Traces.CA(1:end-1),data{4}.Thermo.HRRfilt,'-','COlor', Colour4, 'LIneWidth',1.5);
pl(5) = plot(data{5}.Traces.CA(1:end-1),data{5}.Thermo.HRRfilt,'-','COlor', Colour5, 'LIneWidth',1.5);
pl(6) = plot(data{6}.Traces.CA(1:end-1),data{6}.Thermo.HRRfilt,'-','COlor', Colour6, 'LIneWidth',1.5);
%pl(7) = plot(data{7}.Traces.CA(1:end-1),data{7}.Thermo.HRRfilt,'-','COlor', Colour7, 'LIneWidth',1.5);
%pl(8) = plot(data{8}.Traces.CA(1:end-1),data{8}.Thermo.HRRfilt,'-','COlor', Colour8, 'LIneWidth',1.5);
%pl(9) = plot(data{9}.Traces.CA(1:end-1),data{9}.Thermo.HRRfilt,'-','COlor', Colour9, 'LIneWidth',1.5);
%pl(10) = plot(data{10}.Traces.CA(1:end-1),data{10}.Thermo.HRRfilt,'-','COlor', Colour10, 'LIneWidth',1.5);

%matching the ticks and ylims with the pressure trace plot
xlim([-10 40]);ylim([-12.5 400]);set(gca,'Ytick',0:25:300,'YtickLabel',{'0','','50','','100','','150','','200','','250','','300','','350','',''});
set(gca,'FOntSize',Fsize,'FontWeight','bold');
% set(get(gca,'Xlabel'),'String','CA [°]');
set(get(gca,'Ylabel'),'String','aHRR [J/CA°]');
set(get(gca,'Ylabel'),'Position',get(get(gca,'Ylabel'),'Position')+[0 -50 0]);
%set(get(gca,'Title'),'String','95% H_2 injection with split ratio variation','FontSize',Fsize,'FontWeight','bold');

%Plot the Diesel injection timing markings
TextMarking={'Diesel','10/90/D','20/80/D','30/70/D','40/60/D','50/50/D','60/40/D','70/30/D','80/20/D','90/10/D'};
minSOI=min(SOI);
for i=1:length(data)
    plot(SOI(i)+DOI(i),(82-10*i),'>','Color',get(p(i),'Color'),'MarkerFaceColor',get(p(i),'Color'));
    plot([SOI(i) SOI(i)+DOI(i)],[1 1 ]*(82-10*i),'-','Color',get(p(i),'Color'));
    plot([SOI(i) SOI(i)],[-20,(87-10*i)],'-','Color',get(p(i),'Color'));
    text(minSOI-0.2,82-10*i,TextMarking{i},'HorizontalAlignment','right','FontSize',Fsize-11,'Color',get(p(i),'Color'));
end

%plot the legend. The location of the legend will have to be adjusted
%manually afterwards
lgd = legend(p(1:6),{'Diesel','10/90','20/80','30/70','40/60','50/50','60/40','70/30','80/20','90/10'},'Location','southeast','FontSize',Fsize-9);
lgd.Title.String = sprintf('95%% H_2 injection\nsplit ratio');
lgd.Title.FontSize = Fsize-8;

text(0.13,0.97,'95% H_2 injection with split ratio variation','Units','normalized','FontWeight','bold','FontSize',12)