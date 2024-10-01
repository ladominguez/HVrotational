 clear
format short

%rutaarch = 'C:\Users\mbaen\POSDOC\MBR\RuidoGEOF\';
%rutagrab = 'C:\Users\mbaen\POSDOC\MBR\RuidoGEOFmat\';

% rutaarch = '/Users/antonio/Dropbox/Geofisica/Research/HVrotational/data/';
% rutagrab = '/Users/antonio/Dropbox/Geofisica/Research/HVrotational/matlab/mat/';

cargar_rutas_locales_datos
addpath('utils')

sep = obtener_separador_linux_window();

listest = dir(fullfile(rutaarch));
listest = {listest.name}';
bal = find(ismember(listest,[{'.'};{'..'}])==1);
listest(bal) = [];

rev = [];
cont = 0;
for i = 1:length(listest)
    if isfolder([rutaarch,listest{i}])
        cont = cont+1;
        rev(cont,1) = i;
    end
end
listest = listest(rev);

%% DATOS INICIALES
senhal = 'Ruido';
tiporeg = 'velo'; %acel; velo

% Filtro inicial de las señales en bruto
w1 = 0.01;
w2 = 49.9;

% Factor de esquina para taper de señales
factap = 0.0;

% Elegir nuevo dt
dt = 0.01;

%% S
% buscar = {'NARA'};    %¡¡¡ESCOGER ESTACIÓN!!!%
buscar = listest;

%%
if ~exist(rutagrab,'dir'); mkdir(rutagrab); end

% Ciclo estaciones
[~,Nbuscar] = ismember(buscar,listest);
for k = 2:length(buscar)
    estac = listest{Nbuscar(k)};

    fprintf(1,'%d%s%d%s%s\n',k,sep,length(buscar),' --> ',estac);
    if ~exist([rutagrab,estac,sep,tiporeg],'dir')
        mkdir([rutagrab,estac,sep,tiporeg]);
    end

    listreg = dir(fullfile(rutaarch,estac,'*.sac'));
    listreg = {listreg.name}'; %name
    listregtot = cell(length(listreg),4);
    cont0 = [];
    interv = 100:120;
    for i = interv %length(listreg)
        REG = [rutaarch,estac,sep,listreg{i}];
        [~,~,datREG] = rdsac(REG);
        ymd = datetime(datREG.NZYEAR,01,0)+days(datREG.NZJDAY);
        hh = ['0',num2str(datREG.NZHOUR)];
        mm = ['0',num2str(datREG.NZMIN)];
        fechaymd = [datestr(ymd,'yyyymmdd'),hh(end-1:end),mm(end-1:end)];
        fechajuliana = [num2str(datREG.NZYEAR),num2str(datREG.NZJDAY),hh(end-1:end),mm(end-1:end)];
        if contains(datREG.KCMPNM,'E'); direc = 'EW'; end
        if contains(datREG.KCMPNM,'N'); direc = 'NS'; end
        if contains(datREG.KCMPNM,'Z'); direc = 'VE'; end
        listregtot(i,:) = [listreg(i) {fechajuliana} {fechaymd} {direc}];
        cont0 = [cont0;i];
    end
    listregtot = listregtot(cont0,:);
    dias = unique(listregtot(:,2));

    buscardia = dias;
    % buscardia = {'20232170000';'20232190000'};   % ¡¡¡ESCOGER FECHAS!!!
    [~,Nbuscardia] = ismember(buscardia,dias);
    for ee = 1:length(Nbuscardia)
        m = Nbuscardia(ee);
        fprintf(1,'\t%s%s%d%s%d%s%s\n',estac,' --> ',ee,'/',length(Nbuscardia),': ',dias{m});

        indreg = find(ismember(listregtot(:,2),dias(m))==1);
        fechajuliana = listregtot{indreg(1),2};
        fechaymd = listregtot{indreg(1),3};
        nombarchgrab = [rutagrab,estac,sep,tiporeg,sep,fechaymd,'.mat'];
        % if exist(nombarchgrab,'file') ~= 0; continue; end

        % Lectura SAC
        for ii = 1:3
            if strcmp(listregtot{indreg(ii),4},'EW')
                EWsac = [rutaarch,estac,sep,listregtot{indreg(ii),1}];
                [datEW,iniEW,datosEW] = rdsac(EWsac);
                datEW = double(datEW);
            elseif strcmp(listregtot{indreg(ii),4},'NS')
                NSsac = [rutaarch,estac,sep,listregtot{indreg(ii),1}];
                [datNS,iniNS,datosNS] = rdsac(NSsac);
                datNS = double(datNS);
            elseif strcmp(listregtot{indreg(ii),4},'VE')
                VEsac = [rutaarch,estac,sep,listregtot{indreg(ii),1}];
                [datVE,iniVE,datosVE] = rdsac(VEsac);
                datVE = double(datVE);
            end
        end
        dtorig = datosEW.DELTA;
        NQ = 1/(2*dtorig);

        % Filtro final de las señales en bruto
        w2 = NQ-1;  % 0=s/filtro

        datEW(isnan(datEW)) = 0;
        datNS(isnan(datNS)) = 0;
        datVE(isnan(datVE)) = 0;

        %%%%% EVALUACIÓN %%%%%
        % if length(datEW) < Nseg_sig*mpsVE*0.40 || length(datNS) < Nseg_sig*mpsVE*0.40 || length(datVE) < Nseg_sig*mpsVE*0.40
        %     continue
        % end
        % Nrev = 0.2*length(datEW);
        % if sum(datEW(1+Nrev:end-Nrev)) == 0 || sum(datNS(1+Nrev:end-Nrev)) == 0 || sum(datVE(1+Nrev:end-Nrev)) == 0
        %     continue
        % end
        %%%%%%%%%%%%%%%%%%%%%%

        % % Corrección instrumental
        % [facV,facN,facE] = factorequipo(['T',sensor]);
        % datEW = datEW*facE{1};
        % datNS = datNS*facN{1};
        % datVE = datVE*facV{1};

        % Hora de inicio común para las tres direcciones
        iniDIRvec = datevec(max([iniEW iniNS iniVE]));
        iniEWvec = datevec(iniEW);
        iniNSvec = datevec(iniNS);
        iniVEvec = datevec(iniVE);
        EWini = round(((iniDIRvec(5)*60+iniDIRvec(6))-(iniEWvec(5)*60+iniEWvec(6)))/dt)+1;
        NSini = round(((iniDIRvec(5)*60+iniDIRvec(6))-(iniNSvec(5)*60+iniNSvec(6)))/dt)+1;
        VEini = round(((iniDIRvec(5)*60+iniDIRvec(6))-(iniVEvec(5)*60+iniVEvec(6)))/dt)+1;
        datEW = datEW(EWini:end);
        datNS = datNS(NSini:end);
        datVE = datVE(VEini:end);
        mindat = min([length(datEW) length(datNS) length(datVE)]);
        datEW = datEW(1:mindat);
        datNS = datNS(1:mindat);
        datVE = datVE(1:mindat);

        % Diezmado de la señal
        durac = (length(datEW)-1)*dtorig;
        tnuevo = (0:dt:durac)';
        dtorig = ceil(dtorig*1000)/1000;
        if dtorig ~= dt
            torig = (0:dtorig:durac)';
            datEW = interp1(torig,datEW,tnuevo);
            datNS = interp1(torig,datNS,tnuevo);
            datVE = interp1(torig,datVE,tnuevo);
        end

        % Corta el comienzo y final de la señal (minutos)
        mincutini = dt/60;  % 5
        mincutfin = 0;      % 5
        datEW = datEW(mincutini*60/dt:end-mincutfin*60/dt,:);
        datNS = datNS(mincutini*60/dt:end-mincutfin*60/dt,:);
        datVE = datVE(mincutini*60/dt:end-mincutfin*60/dt,:);
        mindat = min([length(datEW) length(datNS) length(datVE)]);

        % % Figura de revisión
        % figure(20)
        % tnuevo = (0:dt:(length(datEW)-1)*dt)';
        % ml = max([max(abs(datEW)) max(abs(datNS)) max(abs(datVE))]);
        % plot(tnuevo,datNS/ml,'b',tnuevo,datEW/ml+1,'b',tnuevo,datVE/ml+2,'b'); hold on; grid on

        % Quita la tendencia de las señales completas
        datEW = detrend(datEW);
        datNS = detrend(datNS);
        datVE = detrend(datVE);

        % Resta la media de las señales completas
        meandatEW = sum(datEW)/length(datEW);
        meandatNS = sum(datNS)/length(datNS);
        meandatVE = sum(datVE)/length(datVE);
        datEW = datEW-meandatEW;
        datNS = datNS-meandatNS;
        datVE = datVE-meandatVE;

        % Filtro de la señal
        if w1 > 0 && w2 > 0
            datEW = filtsig(datEW,dt,w1,w2,factap);
            datNS = filtsig(datNS,dt,w1,w2,factap);
            datVE = filtsig(datVE,dt,w1,w2,factap);
        end

        % Taper
        EW = datEW.*tukeywin(length(datEW),factap);
        NS = datNS.*tukeywin(length(datNS),factap);
        VE = datVE.*tukeywin(length(datVE),factap);

        % Longitud de los vectores par
        N = length(EW);
        if rem(N,2) ~= 0
            EW = EW(1:end-1);
            NS = NS(1:end-1);
            VE = VE(1:end-1);
        end

        % % Figura de revisión
        % figure(20)
        % t = (0:dt:(length(EW)-1)*dt)';
        % plot(t,NS/ml,'m',t,EW/ml+1,'m',t,VE/ml+2,'m'); hold on; grid on
        % text(t(end)+10,0,'NS'); text(t(end)+10,1,'EW'); text(t(end)+10,2,'VE')
        % set(gcf,'color','white')
        % set(gca,'YTick',[])
        % xlabel('Tiempo (s)')
        % ylabel('Velocidad (cm/s)')
        % ylim([-0.5 2.5])
        % line([t(1) t(end)],[0 0],'color','k')
        % line([t(1) t(end)],[1 1],'color','k')
        % line([t(1) t(end)],[2 2],'color','k')
        % cla

        save(nombarchgrab,'estac','dt','EW','NS','VE','fechajuliana','fechaymd','w1','w2','-v7.3');
        % save(nombarchgrab,'estac','sensor','dt','EW','NS','VE','fechahmsreg','-v7.3');
    end
end
