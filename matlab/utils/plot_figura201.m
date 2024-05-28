function plot_figura201(dt, Smin, Smax, STALTAVE)

    figure(201)
    subplot(2,1,1)
    t =  (0:dt:(length(STALTAVE)-1)*dt).';
    plot(t,STALTAVE,'k'); 
    hold on %; grid on
    line([t(1) t(end)],[Smax Smax],'color','r','linestyle','--','linewidth',2)
    line([t(1) t(end)],[Smin Smin],'color','r','linestyle',':','linewidth',2)
    set(gca,'YTick',0:1:3)
    ylabel('STA/LTA ratio','fontname','Times New Roman','fontSize',14)
    xlim([0 t(end)])
    ylim([0 3.5])
    set(gca,'XTickLabel',[])
    set(gca,'fontname','Times New Roman','fontSize',14)
    h1 = plot(0,0,'k');
    h2 = plot(0,0,'--r','linewidth',2);
    h3 = plot(0,0,':r','linewidth',2);
    lg = legend([h1 h2 h3],'STA/LTA',['(STA/LTA)max = ',num2str(Smax)],['(STA/LTA)min = ',num2str(Smin)]);
    set(lg,'location','northwest','fontname','Times New Roman','fontSize',14)