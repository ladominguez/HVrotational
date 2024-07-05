function [h, leyenda] = figure_ee(ee, f, HV, leyenda, estacion)
    h = figure(ee);
    N = length(HV.clavecomb);
    for ic = 1:N
        % HVmeansmooth = suavmatr(HV.HVmean_comb{ic},HV.fcomb{ic},0,24); %length(HV.HVmean_comb{ic})*0.1
        semilogx(HV.fcomb{ic},HV.HVmean_comb{ic},'linewidth',1.5); hold on; grid on
        % plot(HV.fcomb{ic},HVmeansmooth,'linewidth',1.5); hold on %,'color',col(ic,:) hold on; grid on
        % semilogx(HV.fcomb{ic},HV.NVmean_comb{ic},'linewidth',1.5); hold on; grid on
        % semilogx(HV.fcomb{ic},HV.EVmean_comb{ic},'linewidth',1.5); hold on; grid on
    end
    if isempty(ic); ic = 0; end
    leyenda = [leyenda;HV.clavecomb];
    leg = legend(leyenda);
    title(['HVSR ',estacion],'fontname','Liberation Serif','fontSize',12,'interpreter','none')
    set(gca,'fontname','Liberation Serif','fontSize',12)
    set(gcf,'color','white')
    xlabel('Frecuencia (Hz)','fontname','Liberation Serif','fontSize',12)
    ylabel('Amplitud H/V','fontname','Liberation Serif','fontSize',12)
    maxyplot = get(gca,'ytick');
    xlim([min(f) max(f)]); %[0.01 10]
    % set(gca,'xtick',[0,1:1:10]) %[0.1,1:1:10]