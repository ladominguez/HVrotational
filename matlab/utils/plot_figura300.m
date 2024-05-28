
function plot_figura300(dt, Smax, STALTANS, STALTAEW, STALTAVE)

    figure(300)
    set(gcf,'Position',get(0,'Screensize'))
    fig = tiledlayout(3,2,'TileSpacing','tight','Padding','tight');
    
    t = (0:dt:(length(STALTANS)-1)*dt).';
    
    nexttile(1)
    plot(t,STALTANS,'k'); hold on; grid on
    line([t(1) t(end)],[Smax Smax],'color','r','linestyle','--','linewidth',2)
    set(gca,'YTick',0:1:3)
    set(gca,'XTickLabel',[])
    ylabel('NS','fontname','Times New Roman','fontSize',14)
    ylim([0 11])
    set(gca,'fontname','Times New Roman','fontSize',14)
    
    nexttile(3)
    plot(t,STALTAEW,'k'); hold on; grid on
    line([t(1) t(end)],[Smax Smax],'color','r','linestyle','--','linewidth',2)
    set(gca,'YTick',0:1:3)
    set(gca,'XTickLabel',[])
    ylabel('EW','fontname','Times New Roman','fontSize',14)
    ylim([0 11])
    set(gca,'fontname','Times New Roman','fontSize',14)
    
    nexttile(5)
    plot(t,STALTAVE,'k'); hold on; grid on
    line([t(1) t(end)],[Smax Smax],'color','r','linestyle','--','linewidth',2)
    set(gca,'YTick',0:1:3)
    ylabel('VE','fontname','Times New Roman','fontSize',14)
    xlabel('Time (s)','fontname','Times New Roman','fontSize',14)
    ylim([0 11])
    set(gca,'fontname','Times New Roman','fontSize',14)
    h1 = plot(0,0,'k');
    h2 = plot(0,0,'--r','linewidth',2);
    % lg = legend([h1 h2],'STA/LTA',['(STA/LTA)max = ',num2str(Smax)]);
    % set(lg,'location','south outside','fontname','Times New Roman','fontSize',14)