function [EW,NS,VE,dt,fechahms,w1,w2] = F_dat2mat(EWsac,NSsac,VEsac)
disp(NSsac)
% Lectura SAC
coldat = 2;
[datEW,~,datosEW] = rdsac(EWsac);
[datNS,~,datosNS] = rdsac(NSsac);
[datVE,~,datosVE] = rdsac(VEsac);
datEW = double(datEW);
datNS = double(datNS);
datVE = double(datVE);

ano = datosEW.NZYEAR;
diajuliano = datosEW.NZJDAY;
ymd = datetime(ano,01,0)+days(diajuliano);
fechahms = [datestr(ymd,'yyyymmdd'),'_',num2str(datosEW.NZHOUR),num2str(datosEW.NZMIN)];
dtorig = datosEW.DELTA;
NQ = 1/(2*dtorig);

% Elegir nuevo dt
dt = 0.01;

% Filtro inicial de las señales en bruto
w1 = 0.001;     % 0=s/filtro
w2 = NQ-1;  % 0=s/filtro

% Factor de esquina para taper de señales
factap = 0.0;

datEW(isnan(datEW)) = 0;
datNS(isnan(datNS)) = 0;
datVE(isnan(datVE)) = 0;

% % Corrección instrumental
% [facV,facN,facE] = factorequipo(['T',sensor]);
% datEW = datEW*facE{1};
% datNS = datNS*facN{1};
% datVE = datVE*facV{1};

% Diezmado de la señal
durac = (length(datEW)-1)*dtorig;
tnuevo = (0:dt:durac)';
if dtorig ~= dt
    torig = (0:dtorig:durac)';
    datEW = interp1(torig,datEW,tnuevo);
    datNS = interp1(torig,datNS,tnuevo);
    datVE = interp1(torig,datVE,tnuevo);
end

% Corta el comienzo y final de la señal (minutos)
mincutini = 5; %dt/60;
mincutfin = 5; %0 37
datEW = datEW(mincutini*60/dt:end-mincutfin*60/dt,:);
datNS = datNS(mincutini*60/dt:end-mincutfin*60/dt,:);
datVE = datVE(mincutini*60/dt:end-mincutfin*60/dt,:);

mindat = min([length(datEW) length(datNS) length(datVE)]);
grsig = {1:mindat};
% % Localización de gaps -------------------
% [grsig,grgap] = locgaps(datEW,dt);
% if isempty(grsig); grsig = {1:mindat}; end
% % ----------------------------------------

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
sumdatEW = 0;
sumdatNS = 0;
sumdatVE = 0;
Ngrsig = 0;
for gr = 1:length(grsig)
    ii = grsig{gr}(1);
    ff = grsig{gr}(end);
    sumdatEW = sumdatEW+sum(datEW(ii:ff));
    sumdatNS = sumdatNS+sum(datNS(ii:ff));
    sumdatVE = sumdatVE+sum(datVE(ii:ff));
    Ngrsig = Ngrsig+length(datEW(ii:ff));
end
meandatEW = sumdatEW/Ngrsig;
meandatNS = sumdatNS/Ngrsig;
meandatVE = sumdatVE/Ngrsig;
datEW = datEW-meandatEW;
datNS = datNS-meandatNS;
datVE = datVE-meandatVE;

% Filtro de la señal
if w1 > 0 && w2 > 0
    datEW = filtsig(datEW,dt,w1,w2,factap);
    datNS = filtsig(datNS,dt,w1,w2,factap);
    datVE = filtsig(datVE,dt,w1,w2,factap);
end

datEW2 = zeros(mindat,1);
datNS2 = zeros(mindat,1);
datVE2 = zeros(mindat,1);
if length(grsig) > 1; fprintf(1,'\t%s%s\n','gaps ',date); end
for gr = 1:length(grsig)
    ii = grsig{gr}(1);
    ff = grsig{gr}(end);
    
    % Taper
    factap2 = 0.0;
    datEW2(ii:ff) = datEW(ii:ff).*tukeywin(length(datEW(ii:ff)),factap2);
    datNS2(ii:ff) = datNS(ii:ff).*tukeywin(length(datNS(ii:ff)),factap2);
    datVE2(ii:ff) = datVE(ii:ff).*tukeywin(length(datVE(ii:ff)),factap2);
end

% Longitud de los vectores par
N = length(datEW2);
if rem(N,2) ~= 0
    datEW2 = datEW2(1:end-1); datEW2(end) = 0;
    datNS2 = datNS2(1:end-1); datNS2(end) = 0;
    datVE2 = datVE2(1:end-1); datVE2(end) = 0;
end

EW = datEW2;
NS = datNS2;
VE = datVE2;
% factap2 = 0.0001;
% EW = EW(1:Ndat).*tukeywin(Ndat,factap2);
% NS = NS(1:Ndat).*tukeywin(Ndat,factap2);
% VE = VE(1:Ndat).*tukeywin(Ndat,factap2);

% % Figura de revisión
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
