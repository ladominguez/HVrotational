clear

cargar_rutas_locales
addpath('utils')

listest0 = dir(fullfile(rutagrab));
bal = [listest0.isdir]';
listest0 = {listest0.name}';
listest = listest0(bal);
bal = find(ismember(listest,[{'.'};{'..'}])==1);
listest(bal) = [];

% fid = fopen(rutaestac);
% textscan(fid,'%s',3);
% estacRED = textscan(fid,'%s %f %f %*[^\n]');
% fclose(fid);
% vecest = estacRED{1};
% latest = estacRED{3};
% lonest = estacRED{2};

% buscar = listest;
buscar = {'TOME'};

listabal = cell2mat(listest);
[~,Nbuscar] = ismember(buscar,listabal);
Nest = length(buscar);

figure(100)
hjet = colormap(jet);
close(100)

figure(100)
hhot = colormap(flipud(hot));
hhot = hhot(10:end,:);
close(100)

f0lista = [];
Ntetalista = [];
Ntetalista2 = [];
gamamaxlista = [];
gamaminlista = [];
for k = 1:length(buscar)
    estac = listest{Nbuscar(k)};
    fprintf(1,'%d%s%d%s%s\n',k,'/',length(buscar),' --> ',estac);

    listreg = dir(fullfile(rutagrab,estac,'*.mat'));
    listreg = {listreg.name}';
    listreg = strrep(listreg,'.mat','');

    load(fullfile(rutagrab,estac,listreg{1}));
    tetavec = HV.tetarot;
    ind90 = (length(tetavec)+1)/2;
    HVdir = HV.HVdir_comb1;
    HVdircont = HV.HVdir_comb1([ind90:end,2:ind90],:);
    fgrab = HV.f_comb1;
    HVmeangrab = HV.HVtot_comb1;

    flim1 = fgrab(1); %20
    flim2 = fgrab(end);
    
    Nflim1 = find(fgrab>=flim1,1,'first');
    Nflim2 = find(fgrab>=flim2,1,'first');
    [~,Nfmax0] = max(HVmeangrab(Nflim1:Nflim2));
    Nfmax = Nfmax0+Nflim1-1;
    flim3 = fgrab(Nfmax)-0.0; %0.17
    flim4 = fgrab(Nfmax)+0.0; %0.22
    Nf1 = find(fgrab>=flim3,1); %
    Nf2 = find(fgrab>=flim4,1); %
    flim1graf = 0.0; %20 fgrab(end)  % para graficar
    flim2graf = 50; %55 fgrab(end)  % para graficar
    paso = 1;
    Ntickmax = 6;

    f0lista = [f0lista;fgrab(Nfmax)];
    colfill = 'k';

    Nflim1graf = find(fgrab>=flim1graf,1,'first');
    Nflim2graf = find(fgrab>=flim2graf,1,'first');

    if fgrab(Nf1) < flim1; Nf1 = Nflim1; end
    if fgrab(Nf2) > flim2; Nf2 = Nflim2; end
    gamanumerador = abs((HVdir(:,Nf1:Nf2))-(HVdircont(:,Nf1:Nf2)));
    gamadenominador = zeros(length(tetavec),Nf2-Nf1+1);
    cont = 0;
    for kk = Nf1:Nf2
        cont = cont+1;
        for gg = 1:length(tetavec)
            bal = min([HVdir(gg,kk),HVdircont(gg,kk)]);
            gamadenominador(gg,cont) = bal;
        end
    end
    gamamean = mean(gamanumerador./gamadenominador,2);
    gamamode = mode(gamanumerador./gamadenominador,2);
    gama = gamamean;

    [gamamax0,Ntetamax0] = maximos(gama,2);
    if length(Ntetamax0) == 1
        Ntetamax0 = [Ntetamax0;length(tetavec)];
        gamamax0 = [gamamax0;gama(end)];
    end
    [~,Norden] = sort(gamamax0);
    Norden = Norden(end-1:end);
    gamamax = gama(Ntetamax0(Norden));
    gamamaxlista = [gamamaxlista;gamamax(end)];
%     [~,bal] = max(HVdircont(Ntetamax0(Norden),Nf1:Nf2));
%     Nteta = Ntetamax0(Norden(bal));
    [~,bal1] = max(HVdircont(:,Nf1:Nf2));
    [~,bal2] = max(max(HVdircont(:,Nf1:Nf2)));
    Nteta = bal1(bal2);
    
    [gamamin0,Ntetamin0] = minimos(gama,2);
    if length(Ntetamin0) == 1
        Ntetamin0 = [1;Ntetamin0];
        gamamin0 = [gama(1);gamamin0];
    end
    [~,Norden2] = sort(gamamin0);
    Norden2 = Norden2(1:2);
    gamamin = gama(Ntetamin0(Norden2));
    if abs(tetavec(Ntetamin0(Norden2(1)))-tetavec(Ntetamin0(Norden2(2)))) < 85
        Norden2(2) = [];
        gamamin(2) = [];
        Ntetamin0 = [1;Ntetamin0];
        gamamin0 = [gama(1);gamamin0];
    end
    [~,Norden2] = sort(gamamin0);
    Norden2 = Norden2(1:2);
    gamamin = gama(Ntetamin0(Norden2));
    gamaminlista = [gamaminlista;gamamin(1)];
    Nteta2 = min(Ntetamin0(Norden2));

    Ntetalista = [Ntetalista;Nteta];
    Ntetalista2 = [Ntetalista2;Nteta2];
        
    figure
    tiledlayout(2,2);
    set(gcf,'position',get(0,'Screensize'))
%     if k == 1
%         h = tiledlayout(2,3);%
%     end
    
    nexttile
    HVSRdib = HVdircont;
    contourf(fgrab(Nflim1graf:Nflim2graf),tetavec,HVSRdib(:,Nflim1graf:Nflim2graf),'linecolor','none'); shading interp
    line([flim1graf flim2graf],[90 90],'color','k')
    view([0 0 1])
    % xlim([0.01 fgrab(end)])
    xlim([flim1graf flim2graf])
    ylim([0 180])
    title(estac,'fontname','Times New Roman','fontSize',14,'interpreter','latex')
    xlabel('Frequency, $f$ (Hz)','fontname','Times New Roman','fontsize',14,'interpreter','latex')
    ylabel('Rotation angle, $\phi$ (deg)','fontname','Times New Roman','fontsize',14,'interpreter','latex')
    set(gca,'xscale','log')
    set(gca,'ytick',0:30:180)
    set(gca,'fontname','Times New Roman','fontsize',14)
    set(gcf,'color','white')
    grid on
    colormap(jet)
    cb = colorbar;
    caxis([0 (max(max(HVSRdib(:,Nflim1graf:Nflim2graf))))])
    % cb.Label.String = '$H/V$ amplitude';
    set(cb,'fontname','Times New Roman','fontSize',14)
    text(30,-37,[{'$H^\prime/V(\phi,f)$'};{'amplitude'}],'fontname','Times New Roman', ...
        'fontsize',14,'interpreter','latex','verticalalignment','bottom')

    nexttile
    % for iii = 1:length(tetavec)
    %     Nteta = iii;
    semilogx(fgrab(Nflim1graf:Nflim2graf),HVmeangrab(Nflim1graf:Nflim2graf),'k','linewidth',2); hold on
    semilogx(fgrab(Nflim1graf:Nflim2graf),HVdircont(Nteta,Nflim1graf:Nflim2graf),'r','linewidth',2); hold on
    semilogx(fgrab(Nflim1graf:Nflim2graf),HVdir(Nteta,Nflim1graf:Nflim2graf),'--b','linewidth',2); hold on
    limy = max(HVmeangrab(Nflim1graf:Nflim2graf));
    fill(fgrab([Nf1,Nf2,Nf2,Nf1,Nf1]),[0,0,limy,limy,0],colfill,'facealpha',0.2); hold on
    xlim([flim1graf flim2graf])
    % ylim([0 20]) %[0 limy+3]
    % title([num2str(tetavec(Nteta)),' deg'],'fontname','Times New Roman','fontSize',11)
    xlabel('Frequency, $f$ (Hz)','fontname','Times New Roman','fontsize',14,'interpreter','latex')
    ylabel('$H/V$ amplitude','fontname','Times New Roman','fontsize',14,'interpreter','latex')
    texto = [{['f = ',num2str((round(fgrab(Nfmax)*100)/100)),' Hz']}, ...
        {['γ = ',num2str((round(gamamax(end)*10000)/10000)*100),'%']}];
    text(fgrab(Nfmax),1,texto,'fontname','Times New Roman','fontsize',14)
%     if k == 1 || k == 2
        % text(0.2,1.5,['\gamma = ',num2str((round(gamamax(end)*10000)/10000)*100),'%'],'fontname','Times New Roman','fontsize',14)
%     end
    set(gca,'xscale','log')
    set(gca,'fontname','Times New Roman','fontsize',14)
    set(gcf,'color','white')
    grid on
    tetamod = tetavec(Nteta)+90;
    if tetamod > 180;tetamod = tetamod-180; end
    % leg = legend('$H/V$ total',['$H^\prime/V(\phi,f), \phi=',num2str(tetavec(Nteta)),'^o$'],['$H^\prime/V(\phi,f), \phi=',num2str(tetamod),'^o$'],[num2str(flim3),' Hz$\leq f\leq$',num2str(flim4),' Hz']);
    leg = legend('$H/V$ total',['$H^\prime/V(\phi,f), \phi=',num2str(tetavec(Nteta)),'^o$'],['$H^\prime/V(\phi,f), \phi=',num2str(tetamod),'^o$']);
    set(leg,'fontname','Times New Roman','fontSize',13,'interpreter','latex')
    
    saveas(gcf,[rutagrab,estac,'.png'])
    % print(gcf,[rutagrab,estac],'-dpng','-r300')
        
    listafrec(k,:) = [str2num(estac(3:end)) fgrab(Nfmax)];
end
tetastr = num2str(tetavec(Ntetalista).');

T0lista = 1./f0lista;
maxT0lista = max(T0lista);
minT0lista = min(T0lista);
porc = length(hjet(:,1))/maxT0lista;
