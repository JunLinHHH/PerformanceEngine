    figure;
    
    
    Num_Rows = 2;
    Num_coloum = 3;
    
    To_edge = 65;
    Between_figure_vert = 15;
    Between_figure_hori = 80;

    aaxes_Height = 250;
    aaxes_width = 300;
    
    set(gcf,'unit','pixel','position',[60 60 (To_edge+Between_figure_hori+aaxes_width)+(Num_coloum-1)*360 To_edge+aaxes_Height+Between_figure_vert+(Num_Rows-1)*(aaxes_Height+Between_figure_vert)]);
    

    ax11 = axes('unit','pixels','position',[65 65 aaxes_width aaxes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');
    ax12 = axes('unit','pixels','position',[To_edge+aaxes_width+Between_figure_hori To_edge aaxes_width aaxes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');
    ax13 = axes('unit','pixels','position',[To_edge+2*(aaxes_width+Between_figure_hori) To_edge aaxes_width aaxes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');

    ax21 = axes('unit','pixels','position',[65 To_edge+aaxes_Height+Between_figure_vert aaxes_width aaxes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');
    ax22 = axes('unit','pixels','position',[To_edge+aaxes_width+Between_figure_hori To_edge+aaxes_Height+Between_figure_vert aaxes_width aaxes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');
    ax23 = axes('unit','pixels','position',[To_edge+2*(aaxes_width+Between_figure_hori) To_edge+aaxes_Height+Between_figure_vert aaxes_width aaxes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');
    hold on

    %H₂, CO₂, NOₓ, 1ˢᵗ, 2ⁿᵈ,3ʳᵈ ,4ᵗʰ, [°CA aTDC]

    %plot(a1,Data{1:3,18},Data{1:3,2},'linewidth',1.5,'Marker','o','Markersize',6)
    %plot(a1,Data{1:3,18},Data{4:6,2},'linewidth',1.5,'Marker','o','Markersize',6)
    %plot(a1,Data{1:3,18},Data{7:9,2},'linewidth',1.5,'Marker','o','Markersize',6)
    %plot(a1,Data{1:3,18},Data{10:12,2},'linewidth',1.5,'Marker','o','Markersize',6)
    %plot(a1,Data{1:3,18},Data{13:15,2},'linewidth',1.5,'Marker','o','Markersize',6)

    %plot(ax13,Data{1:3,18},Data{7:9,    16},'linewidth',1.5,'Marker','o','Markersize',6)
    %plot(ax13,Data{1:3,18},Data{10:12,  16},'linewidth',1.5,'Marker','o','Markersize',6)
    %plot(ax13,Data{1:3,18},Data{13:15,  16},'linewidth',1.5,'Marker','o','Markersize',6)
    %plot(ax13,Data{1:3,18},Data{16:18,  16},'linewidth',1.5,'Marker','o','Markersize',6)
    %plot(ax13,Data{1:3,18},Data{19:21,  16},'linewidth',1.5,'Marker','o','Markersize',6)
%
    %plot(ax21,-Data{[4 7 10 13 16 19 22 25 28],17},Data{[4 7 10 13 16 19 22 25 28],    2},'linewidth',1.5,'Marker','o','Markersize',6)

    %plot(gca,Data{1:3,18},Data{22:24,  2},'linewidth',1.5,'Marker','o','Markersize',6)

    %roi = drawpolygon('Parent',gca);
    %P = roi.Position;  
    %ph = patch('Parent',gca, 'XData',P(:,1), 'YData',P(:,2), 'FaceColor','k', 'FaceAlpha',0.18, 'EdgeColor','none','HitTest','off', 'PickableParts','none');


    x=plot(ax23,-Data{[4 31 25 28],17},Data{[4 31 25 28],          16},'LineStyle','-','linewidth',1.5,'Marker','o','Markersize',6,'Color',  [0.00,0.45,0.74])
    y=plot(ax23,-Data{[4 31 25 28],17},Data{[4 31 25 28]+1,        16},'LineStyle','-','linewidth',1.5,'Marker','o','Markersize',6,'Color',  [0.85,0.33,0.10])
    z=plot(ax23,-Data{[4 31 25 28],17},Data{[4 31 25 28]+2,        16},'LineStyle','-','linewidth',1.5,'Marker','o','Markersize',6,'Color',  [0.93,0.69,0.13])
    plot(ax13,-Data{[7 10 13 16 19],17},Data{[7 10 13 16 19],    16},'LineStyle','-','linewidth',1.5,'Marker','o','Markersize',6  ,'Color',[0.00,0.45,0.74])
    plot(ax13,-Data{[7 10 13 16 19],17},Data{[7 10 13 16 19]+1,  16},'LineStyle','-','linewidth',1.5,'Marker','o','Markersize',6,'Color',  [0.85,0.33,0.10])
    plot(ax13,-Data{[7 10 13 16 19],17},Data{[7 10 13 16 19]+2,  16},'LineStyle','-','linewidth',1.5,'Marker','o','Markersize',6,'Color',  [0.93,0.69,0.13])

    bar(ax23,-Data{[4 31 25 28],17},Data{[4 31 25 28],          5},'FaceColor',  [0.00,0.45,0.74],BarWidth=0.25)
    bar(ax23,-Data{[4 31 25 28],17},Data{[4 31 25 28]+1,        5},'FaceColor',  [0.85,0.33,0.10],BarWidth=0.25)
    bar(ax23,-Data{[4 31 25 28],17},Data{[4 31 25 28]+2,        5},'FaceColor',  [0.93,0.69,0.13],BarWidth=0.25)
    
    bar(ax13,-Data{[7 10 13 16 19],17},Data{[7 10 13 16 19],    5},'FaceColor',  [0.00,0.45,0.74],BarWidth=0.25)
    bar(ax13,-Data{[7 10 13 16 19],17},Data{[7 10 13 16 19]+1,  5},'FaceColor',  [0.85,0.33,0.10],BarWidth=0.25)
    bar(ax13,-Data{[7 10 13 16 19],17},Data{[7 10 13 16 19]+2,  5},'FaceColor',  [0.93,0.69,0.13],BarWidth=0.25)

    plot(gca,xlim,[Data{1,5} Data{1,5}],'LineStyle','--','linewidth',1.5,'Marker','none' ,'Color',[0.00,0.45,0.74])
    plot(gca,xlim,[Data{2,5} Data{2,5}],'LineStyle','--','linewidth',1.5,'Marker','none' ,'Color',[0.85,0.33,0.10])
    plot(gca,xlim,[Data{3,5} Data{3,5}],'LineStyle','--','linewidth',1.5,'Marker','none' ,'Color',[0.93,0.69,0.13])