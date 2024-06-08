function [ESTR, teta] = rotar_sismogramas(ESTR, teta, Ndias)


for p = 1:Ndias
    ESTR.EWrot{p} = ESTR.EW{p}*cosd(teta)+ESTR.NS{p}*sind(teta);  %longitudinal
    ESTR.NSrot{p} =-ESTR.EW{p}*sind(teta)+ESTR.NS{p}*cosd(teta);  %transversal
end

%% codigo anterior -- borrar
%teta = tetarot(Nteta);
%EW = ESTR.EW*cosd(teta)+ESTR.NS*sind(teta);  %longitudinal
%NS =-ESTR.EW*sind(teta)+ESTR.NS*cosd(teta);  %transversal
%VE = ESTR.VE;
