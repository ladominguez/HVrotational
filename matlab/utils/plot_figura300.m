function plot_figura300(NS,EW,VE,dt,wincleantot,iv,fv, Smax, STALTANS, STALTAEW, STALTAVE)

figure(300)
set(gcf,'Position',get(0,'Screensize'))
fig = tiledlayout(3,2,'TileSpacing','tight','Padding','tight');

t = (0:dt:(length(NS)-1)*dt).';

% Figuras con los cocientes STA/LTA de cada dirección
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

% Figuras de la señal en cada dirección indicando las ventanas efectivas
d = find(wincleantot~=0);
ml = max([max(abs(NS)) max(abs(EW)) max(abs(VE))]);
% ml = 1;

nexttile(2)
plot(t,NS/ml,'k'); hold on; grid on
ylabel('NS','fontname','Times New Roman','fontSize',14)
xlim([t(1) t(end)])
ylim([-1 1])
% set(gca,'YTick',[-0.02,0,0.02])
set(gca,'XTickLabel',[])
set(gca,'fontname','Times New Roman','fontSize',14)

nexttile(4)
plot(t,EW/ml,'k'); hold on; grid on
ylabel('EW','fontname','Times New Roman','fontSize',14)
xlim([t(1) t(end)])
ylim([-1 1])
% set(gca,'YTick',[-0.02,0,0.02])
set(gca,'XTickLabel',[])
set(gca,'fontname','Times New Roman','fontSize',14)

nexttile(6)
plot(t,VE/ml,'k'); hold on; grid on
ylabel('VE','fontname','Times New Roman','fontSize',14)
xlabel('Time (s)','fontname','Times New Roman','fontSize',14)
xlim([t(1) t(end)])
ylim([-1 1])
% set(gca,'YTick',[-0.02,0,0.02])
set(gca,'fontname','Times New Roman','fontSize',14)

for q = 1:length(wincleantot)

    if ismember(q,d)
        NSd = NS(iv(q):fv(q));
        EWd = EW(iv(q):fv(q));
        VEd = VE(iv(q):fv(q));

        nexttile(2)
        fill([t(iv(q)),t(fv(q)),t(fv(q)),t(iv(q)),t(iv(q))], ...
            [-1,-1,1,1,-1],'b','edgecolor','b','facealpha',0.5) ; hold on

        nexttile(4)
        fill([t(iv(q)),t(fv(q)),t(fv(q)),t(iv(q)),t(iv(q))], ...
            [-1,-1,1,1,-1],'b','edgecolor','b','facealpha',0.5) ; hold on

        nexttile(6)
        fill([t(iv(q)),t(fv(q)),t(fv(q)),t(iv(q)),t(iv(q))], ...
            [-1,-1,1,1,-1],'b','edgecolor','b','facealpha',0.5) ; hold on
        % plot(t(iv(q):fv(q)),NSd/ml,'c',t(iv(q):fv(q)), ...
        %     EWd/ml+1,'c',t(iv(q):fv(q)),VEd/ml+2,'c'); hold on; grid on
        % text(t(iv(q)),0.03,num2str(q),'fontname','Times New Roman','fontSize',14)
        set(gcf,'color','white')
        % h1 = plot(0,0,'k');
        % h2 = fill([0,0,0,0,0],[0,0,0,0,0],'b','edgecolor','r','facealpha',0.1);
        % lg = legend([h1 h2],'Corrected signals','Windows selected');
        % set(lg,'location','south outside','fontname','Times New Roman','fontSize',14)
    end
end
