    %% Creating axes for plotting (pressure&aHRR, combustion phasing&burn duration, performance data(exchlude efficiency&Peak P))

    figure(1);
    set(figure(1),'unit','pixel','position',[60 60 1060 360]);
    
    % Each Axes width and height
    Num_Rows = 1;
    Num_coloum = 5;
    axes_Height = 400;
    axes_width = 280;
    set(figure(1),'unit','pixel','position',[60 60 65*2+30*(Num_coloum-1)+Num_coloum*axes_width 65*2+30+Num_Rows*axes_Height]);

    
    % First line
    ax11 = axes('unit','pixels','position',[65 65 axes_width axes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');
    ax12 = axes('unit','pixels','position',[axes_width+65+30 65 axes_width axes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');
    ax13 = axes('unit','pixels','position',[axes_width+65+30+axes_width+30 65 axes_width axes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');
    ax14 = axes('unit','pixels','position',[axes_width+65+30+axes_width+30+axes_width+30 65 axes_width axes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');
    ax15 = axes('unit','pixels','position',[axes_width+65+30+axes_width+30+axes_width+30+axes_width+30 65 axes_width axes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');



    % Second line
    %ax21 = axes('unit','pixels','position',[65 395 axes_width axes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');
    %ax22 = axes('unit','pixels','position',[305 395 axes_width axes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');
    %ax23 = axes('unit','pixels','position',[545 395 axes_width axes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');
    %ax24 = axes('unit','pixels','position',[785 395 axes_width axes_Height],'box','on', 'linewidth',1.2, 'fontsize',12,'fontweight','bold');

    


    
    
   




