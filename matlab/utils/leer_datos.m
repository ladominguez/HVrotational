
%function [ESTR, vecfechahms] = leer_datos(rutaarch, estac, unidad, listreg, listdias)
function [ESTR, w1, w2] = leer_datos(rutaarch, estac, unidad, listreg, NdiasHV, inddia, w1new, w2new)
        sep = obtener_separador_linux_window();
        ESTR = [];
        for i = 1:NdiasHV
            Nreg = inddia+(i-1);
            REG = load([rutaarch,estac,sep,unidad,sep,listreg{Nreg}]);
            ESTR.EW{i} = double(REG.EW);
            ESTR.NS{i} = double(REG.NS);
            ESTR.VE{i} = double(REG.VE);
            ESTR.dt(i) = REG.dt;
            ESTR.vecfechahms{i} = [listreg{Nreg}(1:8),'_',listreg{Nreg}(9:end-4)];
            if i == 1
                dt = REG.dt;
                w1 = REG.w1;
                w2 = REG.w2;
            end

            % Nuevo filtro entre w1 y w2
            if w1new > 0 || w2new > 0
                if w1new ~= 0; w1 = w1new; end
                if w2new ~= 0; w2 = w2new; end
                ESTR.EW{i} = filtsig(ESTR.EW{i},dt,w1,w2,0.0);
                ESTR.NS{i} = filtsig(ESTR.NS{i},dt,w1,w2,0.0);
                ESTR.VE{i} = filtsig(ESTR.VE{i},dt,w1,w2,0.0);
            end
        end


%% <<CODIGO ANTERIOR>> - borrar en el futuro

%cont0 = 0;
%ESTR = [];
%vecfechahms = [];
%separador = obtener_separador_linux_window();

%for i = 1:length(listreg)
%    if strcmp(listreg{i}(1:8),listdias)
%        cont0 = cont0+1;
%        REG = load([rutaarch,estac,separador,unidad,separador,listreg{i}]);
%        ESTR.EW(:,cont0) = double(REG.EW);
%        ESTR.NS(:,cont0) = double(REG.NS);
%        ESTR.VE(:,cont0) = double(REG.VE);
%        ESTR.dt(cont0,1) = REG.dt;
%        vecfechahms = [vecfechahms;{[listreg{i}(1:8),'_',listreg{i}(9:end-4)]}];
%    end
%end

