
function [ESTR, vecfechahms] = leer_datos(rutaarch, estac, unidad, listreg, listdias)

cont0 = 0;
ESTR = [];
vecfechahms = [];
separador = obtener_separador_linux_window();

for i = 1:length(listreg)
    if strcmp(listreg{i}(1:8),listdias)
        cont0 = cont0+1;
        REG = load([rutaarch,estac,separador,unidad,separador,listreg{i}]);
        ESTR.EW(:,cont0) = double(REG.EW);
        ESTR.NS(:,cont0) = double(REG.NS);
        ESTR.VE(:,cont0) = double(REG.VE);
        ESTR.dt(cont0,1) = REG.dt;
        vecfechahms = [vecfechahms;{[listreg{i}(1:8),'_',listreg{i}(9:end-4)]}];
    end
end

