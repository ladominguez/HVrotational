clear; clc
format short

cargar_rutas_locales
% addpath('utils')

%% DATOS INICIALES
senhal = 'noise';
unidad = 'velo';  % DESP VELO ACEL

% Filtro inicial de las se�ales
w1new = 0; %1;     0=s/filtro
w2new = 0; %255.9; 0=s/filtro

% Factor de esquina para taper de se�ales
factap = 0.01;

% *****NORMALIZACI�N*****
% band:   0=ninguna, 2=suma3direcc, 3=SW
% onebit: 1=SI, 0=NO

% SELECCIONAR DATOS
NdiasHV = 3;
segvent = [500];         % Segundos de las ventanas para inversi�n
porctrasl = [25];        % Porcentaje de traslape de las ventanas
normalizac = [2 0];      % Normalizaci�n: [band,onebit]
tiempoHV = [48*60];      % Tiempo (minutos) para c�lculo de cada H/V (Tiempo de registro manipulable)
ventaleatHV = 0;         % 1=ventanas aleatoria, 0=ventanas continuas
NvBootstrap = 1;         % N�mero de ventanas para el boostrap
tSTA = 1; %1.35;         % En segundos
tLTA = 60;               % En segundos
Smax = 5;                % 0=todas las ventanas
Smin = 0.2;
dfnew = 1;

% Si baja el tLTA es m�s conservador
itertot = length(segvent)*length(porctrasl)*length(normalizac(:,1))*length(tiempoHV);
separador = obtener_separador_linux_window();

%% Buscar estaci�n
listest = dir(fullfile(rutaarch));
listest = {listest.name}';
bal = find(ismember(listest,[{'.'};{'..'}])==1);
listest(bal) = [];

buscar = listest;
buscar = {'TOME'};        % ���ESCOGER ESTACI�N!!!

%% Invierte la escala de colores,se puede comentar
col = get_colors(itertot);

%% Ciclo principal
% tetarot = 0:45:180;
tetarot = 0;
[~,Nbuscar] = ismember(buscar,listest);

for aleat = 1 %:10   % se removio este ciclo en la siguiente versi�n

    % Ciclo para estaciones
    for ee = 1:length(buscar)
        estac = listest{Nbuscar(ee)};
        fprintf(1,'%d%s%d%s%s\n',ee,'/',length(buscar),' --> ',estac);

        crear_directorios(rutagrab, estac)
        nombgrab = [rutagrab,estac,[separador 'HV_'],estac];


        listreg = dir(fullfile(rutaarch,estac,unidad,'*.mat'));
        listreg = {listreg.name}'; %name


        listdias = obtener_lista_dias(listreg);

        % buscardia = listdias;
        % buscardia = {'20200929';'20200116';'20201129'};
        % [~,Nbuscardia] = ismember(buscardia,listdias);

        diaini = (1:NdiasHV:length(listdias)).';

        leyenda = [];
        % Loop por cada d�a
        %for dd = 1:length(Nbuscardia)
        for dd = 1:length(diaini)

            inddia = diaini(dd);
            %k = Nbuscardia(dd);
            nombgrab0 = [nombgrab,'_',listdias{inddia},'.mat'];
            % if exist(nombgrab0,'file') ~= 0; continue; end

            %[ESTR, vecfechahms] = leer_datos(rutaarch, estac, unidad, listreg, listdias{k});
            [ESTR, dt, w1, w2] = leer_datos(rutaarch, estac, unidad, listreg, NdiasHV, inddia, w1new, w2new);
            
            dt = ESTR.dt(1);
            fmax = 1/(2*dt);
            
            [~,Ndias] = size(ESTR.EW);
            f0 = [];

            % estac: nombre de la estaci�n
            % paraadic: par�metros adicionales (fechas, n�mero de ventanas para H/V, par�metros para STA/LTA, df)
            % clavecomb: clave de cada combinaci�n de par�metros
            % Nvent: n�mero de ventanas empleadas para el H/V
            % fcomb: vector de frecuencias
            % HVmean_comb: matriz con el H/V medio de cada combinaci�n de par�metros (por columna)
            % NVmean_comb: matriz con el H/V de cada combinaci�n de par�metros, usando solo la componente norte-sur (por columna)
            % EVmean_comb: matriz con el H/V de cada combinaci�n de par�metros, usando solo la componente este-oeste (por columna)
            % tiempoHV_orig_min: tiempo solicitado para el c�lculo del H/V
            % tiempoHV_real_min: tiempo real empleado para el c�lculo del H/V
            % f_comb1: vector de frecuencias de la combinaci�n de par�metros 1
            % HVtot_comb1: H/V de la combinaci�n de par�metros 1
            % HVNSdir_comb1: H/V direccional norte-sur empleando la combinaci�n de par�metros 1
            % HVEWdir_comb1: H/V direccional este-oeste empleando la combinaci�n de par�metros 1
            % tetarot: vector de �ngulos de rotaci�n para el H/V direccional

            HV = struct('estac',[],'paraadic',[],'clavecomb',[],'Nvent',[],'fcomb',[],'HVmean_comb',[],'NVmean_comb',[],'EVmean_comb',[], ...
                'tiempoHV_orig_min',[],'tiempoHV_real_min',[], ...
                'f_comb1',[],'HVtot_comb1',[],'HVNSdir_comb1',[],'HVEWdir_comb1',[],'tetarot',[]);

            HV.HVNSdir_comb1 = [];
            HV.HVEWdir_comb1 = [];

            for Nteta = 1:length(tetarot)
                iter = 0;
                ccd = 0;
                Nv = 0;

                % Rotaci�n sismogramas
                [ESTR teta] = rotar_sismogramas(ESTR, tetarot( Nteta), Ndias );

                fprintf(1,'\t%d%s%d%s%d\n',Nteta,'/',length(tetarot),' --> teta=',teta);

                % *****************************************************
                % CICLO LONGITUD DE VENTANAS
                % *****************************************************
                for vv = 1:length(segvent)

                    [f, fin, ini, ptosvent, Nespec, df] = obtener_vector_de_frecuencia(segvent(vv), ...
                        dt, dfnew, fmax);

                    Nfrecred = fin-ini+1;

                    % *****************************************************
                    % CICLO TRASLAPE DE VENTANAS
                    % *****************************************************
                    wincleantot = [];
                    for tt = 1:length(porctrasl)

                        % VENTANEO
                       
                        [Nventefec, M, iv, fv, wincleantot, wincleanEW, wincleanNS, ...
                            wincleanVE, STALTAEW, STALTANS, STALTAVE] = ventaneo(porctrasl(tt), ...
                            ptosvent, ESTR, dt, tSTA, tLTA, Smax, Smin, Ndias );

                        % % Figuras para revisi�n 1
                        % plot_figura300(NS,EW,VE,dt,wincleantot,iv,fv, Smax, STALTANS, STALTAEW, STALTAVE)
                        % close(300)
                        % plot_figura201(dt, Smin, Smax, STALTAVE)

                        % DIVISI�N DE LA SE�AL EN VENTANAS DE TIEMPO

                        [EWv, NSv, VEv, fechahmsvent ] = division_ventanas_tiempo(ESTR, ptosvent, ...
                            Nventefec, Ndias, wincleantot, iv, fv);

                        % % VENTANAS EFECTIVAS
                        % [EWv, NSv, VEv, fechahmsvent] = ventanas_efectivas(ESTR, ptosvent, Nventefec, Ndias, wincleantot, iv, fv);


                        % REMUEVE LA MEDIA POR VENTANAS Y APLICA TAPER
                        [EWv, NSv, VEv] = remover_media_taper(EWv, NSv, VEv, ptosvent, factap);

                        % % Figuras para revisi�n 2
                        % figure(201)
                        % t = (0:dt:(length(NS)-1)*dt).';
                        % d = find(wincleanVE~=0);
                        % ml = max([max(abs(NS)) max(abs(EW)) max(abs(VE))]);
                        % 
                        % subplot(2,1,2)
                        % plot(t,VEm,'k'); hold on %; grid on
                        % ylabel('Velocity (cm/s)','fontname','Times New Roman','fontSize',14)
                        % xlabel('Time (s)','fontname','Times New Roman','fontSize',14)
                        % ylim([-max(abs(VEm)) max(abs(VEm))])
                        % set(gca,'fontname','Times New Roman','fontSize',14)
                        % for q = d.'
                        %     if q == d(1)
                        %         yampl = 5*max([max(abs(NSm(iv(q):fv(q)))) ...
                        %             max(abs(EWm(iv(q):fv(q)))) max(abs(VEm(iv(q):fv(q))))]);
                        %     end
                        %     fill([t(iv(q)),t(fv(q)),t(fv(q)),t(iv(q)),t(iv(q))], ...
                        %         [-yampl,-yampl,yampl,yampl,-yampl],'b','edgecolor','r','facealpha',0.3) ; hold on
                        % end
                        % xlim([0 t(end)])
                        % ylim([-yampl yampl])
                        % set(gcf,'color','white')
                        % h1 = plot(0,0,'k');
                        % h2 = fill([0,0,0,0,0],[0,0,0,0,0],'b','edgecolor','r','facealpha',0.1);
                        % lg = legend([h1 h2],'Corrected signals','Windows selected');
                        % set(lg,'location','northwest','fontname','Times New Roman','fontSize',14)
                        % % cla

                        % *****************************************************
                        % CICLO DE NORMALIZACI�N
                        % *****************************************************
                        for norm = 1:length(normalizac(:,1))
                            band = normalizac(norm,1);
                            onebit = normalizac(norm,2);

                            % Carpeta y archivo para grabar resultados

                            nombcomb = nombre_combinac(senhal,unidad,band,w1,w2,onebit,segvent(vv),porctrasl(tt));

                            % NORMALIZACI�N
                            [fNSventnorm,fVEventnorm,fEWventnorm,~,~,~,~,~] = ...
                                F_normalizacionfrec(NSv,VEv,EWv, ...
                                Nespec,band,onebit,dt,factap);

                            % promNS_EW = abs((fNSventnorm+fEWventnorm)/2);
                            fNSvent = abs(fNSventnorm(ini:fin,:));
                            fEWvent = abs(fEWventnorm(ini:fin,:));
                            fVEvent = abs(fVEventnorm(ini:fin,:));
                            fHHvent = abs(sqrt((fNSventnorm(ini:fin,:).^2+fEWventnorm(ini:fin,:).^2)/2));

                            fNSventnorm = [];
                            fEWventnorm = [];
                            fVEventnorm = [];

                            % *****************************************************
                            % CICLO DE TIEMPOS PARA C�LCULO DE H/V
                            % *****************************************************
                            for nh = 1:length(tiempoHV)
                                iter = iter+1;
                                fprintf(1,'\t%s%s%s%s%d%s%d\n',estac,'_',listdias{inddia},' --> iter ',iter,'/',itertot);
                                suav = 0;   %0=no; 1=s�

                                [HVtot,NVmean,EVmean,NventHV,vini, ...
                                    tiempoHVnuevo,numHV,HVvent] = ...
                                    F_HVruido(f,fNSvent,fEWvent,fVEvent, ...
                                    fHHvent,segvent(vv),porctrasl(tt), ...
                                    tiempoHV(nh),suav,ventaleatHV,NvBootstrap);

                                tiempoHVnuevo_str = num2str(round(tiempoHVnuevo*100)/100);
                                clavecomb = ['CD-HV',tiempoHVnuevo_str,'hr','-',nombcomb,'-Nw',num2str(NventHV(1)),'-NwBS',num2str(numHV)];

                                if teta == 0
                                    ccd = ccd+1;
                                    HV.clavecomb{ccd} = clavecomb;
                                    HV.Nvent{ccd} = NventHV(1);
                                    HV.fcomb{ccd} = f;
                                    HV.HVmean_comb{ccd} = HVtot;
                                    HV.NVmean_comb{ccd} = NVmean;
                                    HV.EVmean_comb{ccd} = EVmean;
                                    HV.tiempoHV_orig_min{ccd} = tiempoHV(nh);
                                    HV.tiempoHV_real_min{ccd} = tiempoHVnuevo*60;
                                end

                                if teta == 0 && iter == 1
                                    HV.estac = estac;
                                    HV.paraadic.fechahms = ESTR.vecfechahms;
                                    HV.paraadic.ventaleatHV = ventaleatHV;
                                    HV.paraadic.NvBootstrap = NvBootstrap;
                                    HV.paraadic.tSTA = tSTA;
                                    HV.paraadic.tLTA = tLTA;
                                    HV.paraadic.Smax = Smax;
                                    HV.paraadic.Smin = Smin;
                                    if df > dfnew
                                        HV.paraadic.df = dfnew;
                                    else
                                        HV.paraadic.df = df;
                                    end
                                    HV.f_comb1 = HV.fcomb{1};
                                    HV.HVtot_comb1 = HV.HVmean_comb{1};
                                    HV.tetarot = tetarot;
                                end

                                if iter == 1
                                    HV.HVNSdir_comb1 = [HV.HVNSdir_comb1;NVmean.'];
                                    HV.HVEWdir_comb1 = [HV.HVEWdir_comb1;EVmean.'];
                                end
                            end
                            NSv = [];
                            VEv = [];
                            EWv = [];
                        end
                    end
                end

                %% Figura
                if teta == 0
                    h = figure(201);
                    for ic = 1:length(HV.clavecomb)
                        % HVmeansmooth = suavmatr(HV.HVmean_comb{ic},HV.fcomb{ic},0,24); %length(HV.HVmean_comb{ic})*0.1
                        semilogx(HV.fcomb{ic},HV.HVmean_comb{ic},'linewidth',1.5); hold on; grid on
                        % plot(HV.fcomb{ic},HVmeansmooth,'linewidth',1.5); hold on %,'color',col(ic,:) hold on; grid on
                        % semilogx(HV.fcomb{ic},HV.NVmean_comb{ic},'linewidth',1.5); hold on; grid on
                        % semilogx(HV.fcomb{ic},HV.EVmean_comb{ic},'linewidth',1.5); hold on; grid on
                    end
                    if isempty(ic); ic = 0; end
                    leyenda = [leyenda;HV.clavecomb];
                    leg = legend(leyenda);
                    title(['HVSR ',estac],'fontname','Liberation Serif','fontSize',12,'interpreter','none')
                    set(gca,'fontname','Liberation Serif','fontSize',12)
                    set(gcf,'color','white')
                    xlabel('Frecuencia (Hz)','fontname','Liberation Serif','fontSize',12)
                    ylabel('Amplitud H/V','fontname','Liberation Serif','fontSize',12)
                    maxyplot = get(gca,'ytick');
                    xlim([min(f) max(f)]); %[0.01 10]
                    % set(gca,'xtick',[0,1:1:10]) %[0.1,1:1:10]
                    drawnow

                    % print(gcf,nombgrab0(1:end-4),'-dpng','-r600')
                    % close(h)
                end
            end
            save(nombgrab0,'HV','-v7.3');

            % % Figura HV direccional
            % contourf(HV.fcomb{1},HV.tetarot,HV.HVNSdir_comb1); shading interp
            % xlabel('Frecuencia (Hz)','fontname','Liberation Serif','fontSize',12)
            % ylabel('�ngulo de rotaci�n (deg)','fontname','Liberation Serif','fontSize',12)
            % view([0 0 1])
            % xlim([0.01 HV.fcomb{1}(end)])
            % ylim([0 180])
            % colormap(jet)
            % 
            % legtetaevec = num2str(HV.tetarot.');
            % legend(legtetaevec,'interpreter','none')
        end
    end
end
