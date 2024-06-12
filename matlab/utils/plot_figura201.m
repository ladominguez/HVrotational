function plot_figura201(ESTR, dt, Smin, Smax, winclean, STALTA, iv, fv)

for p = 1:length(ESTR.VE)
    figure(201)
    set(gcf,'Position',get(0,'Screensize'))
    tiledlayout(2,1);

    t =  (0:dt:(length(ESTR.VE{p})-1)*dt).';

    nexttile(1)
    plot(t,STALTA.VE{p},'k'); hold on
    line([t(1) t(end)],[Smax Smax],'color','r','linestyle','--','linewidth',2)
    line([t(1) t(end)],[Smin Smin],'color','r','linestyle',':','linewidth',2)
    % set(gca,'YTick',0:1:3)
    ylabel('STA/LTA ratio','fontname','Arial','fontSize',20)
    xlim([t(1) t(end)])
    ylim([0 Smax*2])
    set(gca,'XTickLabel',[])
    set(gca,'fontname','Arial','fontSize',14)
    h1 = plot(0,0,'k');
    h2 = plot(0,0,'--r','linewidth',2);
    h3 = plot(0,0,':r','linewidth',2);
    lg = legend([h1 h2 h3],'STA/LTA',['(STA/LTA)max = ',num2str(Smax)],['(STA/LTA)min = ',num2str(Smin)]);
    set(lg,'location','northwest','fontname','Times New Roman','fontSize',14)

    VEm = ESTR.VE{p};
    % ml = max(max(abs(VEm)));
    ml = 1;
    limy = max(mean(abs(VEm)))*10;
    d = find(winclean.VE{p}~=0);

    nexttile(2)
    plot(t,VEm,'k'); hold on %; grid on
    ylabel('Velocity (cm/s)','fontname','Times New Roman','fontSize',14)
    xlabel('Time (s)','fontname','Times New Roman','fontSize',14)
    xlim([t(1) t(end)])
    ylim([-limy limy])
    set(gca,'fontname','Times New Roman','fontSize',14)
    for q = d.'
        if q == d(1)
            yampl = 5*max(abs(VEm(iv{p}(q):fv{p}(q))));
        end
        fill([t(iv{p}(q)),t(fv{p}(q)),t(fv{p}(q)),t(iv{p}(q)),t(iv{p}(q))], ...
            [-yampl,-yampl,yampl,yampl,-yampl],'b','edgecolor','r','facealpha',0.3) ; hold on
    end
    xlim([0 t(end)])
    ylim([-yampl yampl])
    set(gcf,'color','white')
    h1 = plot(0,0,'k');
    h2 = fill([0,0,0,0,0],[0,0,0,0,0],'b','edgecolor','r','facealpha',0.1);
    lg = legend([h1 h2],'Corrected signals','Windows selected');
    set(lg,'location','northwest','fontname','Times New Roman','fontSize',14)

    pause()
    close(201)
end
