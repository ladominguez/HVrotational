% Elimina picos en la señal
function [winclean,STALTA] = picossig6(SIG,dt,vi,vf,tSTA,tLTA,Smax,Smin)
% SIG = ESTR.EW_mat(:,p);
% dt = 0.01;
% Smax = 3;

SIGabs = abs(SIG);
% SIGabs = abs(hilbert(SIG));
M = length(vi);
N = length(SIG);

winclean = ones(M,1);

% Máximos y medias por ventana
LTAmax = zeros(M,1);
LTAmean = zeros(M,1);
for i = 1:M
    SIGp = SIGabs(vi(i):vf(i));
    LTAmax(i,1) = max(SIGp);
    LTAmean(i,1) = mean(SIGp);
end
LTAmaxmax = max(SIGabs);

% Elimina picos superiores al 99.5% del pico máximo de la señal
for i = 1:M
    if LTAmax(i) > 0.995*LTAmaxmax
        winclean(i) = 0;
    end
end
LTAmeanmax = max(LTAmean(and(LTAmean>0,winclean==1)));
if isempty(LTAmeanmax)
    winclean(winclean==1) = 0;
end
% Elimina picos superiores al 90% del promedio máximo de la señal
for i = 1:M
    if winclean(i) == 1
        if LTAmean(i) > 0.9*LTAmeanmax
            winclean(i) = 0;
        end
    end
end

% Calcula STA y LTA en toda la señal
STA = tSTA/dt;
LTA = tLTA/dt;
STAmovmean = movmean(SIGabs,STA,'omitnan'); %,'Endpoints','discard'
LTAmovmean = movmean(SIGabs,LTA,'omitnan'); %,'Endpoints','discard'
% STAmovmean = [movmean(SIGabs(1:STA-1),STA-1,'omitnan');movmean(SIGabs,STA,'omitnan','Endpoints','discard')];
% LTAmovmean = [movmean(SIGabs(1:LTA-1),LTA-1,'omitnan');movmean(SIGabs,LTA,'omitnan','Endpoints','discard')];
STALTA = STAmovmean./LTAmovmean;

% Evaluación LTAmean/STA > Smax
for i = 1:M
    if winclean(i) == 1
        if max(STALTA(vi(i):vf(i))) >= Smax || max(STALTA(vi(i):vf(i))) <= Smin
            winclean(i) = 0;
        end
    end
end

% Evaluación LTAmean/STA no sea NaN
for i = 1:M
    if winclean(i) == 1
        if isnan(STALTA(vi(i):vf(i)))
            winclean(i) = 0;
        end
    end
end

% figure
% t = (0:dt:(N-1)*dt).';
% subplot(2,1,1)
% plot(t,STALTA)
% line([t(1) t(end)],[Smax Smax],'color','g')
% 
% subplot(2,1,2)
% plot(t,SIG,'b'); hold on
% for i = 1:M
%     if winclean(i) == 1
%         tp = t(vi(i):vf(i));
%         SIGp = SIG(vi(i):vf(i));
%         plot(tp,SIGp,'r'); hold on
%     end
% end
% plot(t,STAmovmean,'k',t,LTAmovmean,'y','linewidth',1.5); hold on
% 1;
