function [nombcomb] = nombre_combinac(senhal,unidad,band,w1,w2,onebit,segvent,porctrasl)

nombcomb = [senhal,'-',unidad,'-','filt',num2str(w1),'to',num2str(w2), ...
    'Hz-band',num2str(band),'-onebit',num2str(onebit),'-tv', ...
    num2str(round(segvent*10000)/10000),'s-Overlap',num2str(porctrasl),'%'];
nombcomb = strrep(nombcomb,' ','');
nombcomb = strrep(nombcomb,'-onebit0','');
nombcomb = strrep(nombcomb,'onebit1','1bit');
nombcomb = strrep(nombcomb,'band0','raw');
nombcomb = strrep(nombcomb,'band1','Eunit');
nombcomb = strrep(nombcomb,'band2','SW3dir');
nombcomb = strrep(nombcomb,'band3','SW');

% % REVISAR BLOQUE ---------
% carpeta = [rutagrab,estac];
% % if ~exist(carpeta,'dir'); mkdir(carpeta); end
% nombarch = [carpeta,'\',estac,'-',nombcomb,'.mat'];
% nombarch = strrep(nombarch,'-','_');
% nombarch = strrep(nombarch,'%','trasl');
% % if exist(nombarch,'file') ~= 0; continue; end
% % ----------------------------------------
