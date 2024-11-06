function [aa,bb] = maxmin(SIG,band)
% band = 1 --> maximos
% band = -1 --> minimos

SIGband = abs(SIG);

if band == 1
    cont = 0;
    for i = 2:length(SIGband)-1
        if SIGband(i) >= SIGband(i-1) && SIGband(i) > SIGband(i+1)
            if SIGband(i) > 0
                cont = cont+1;
                aa(cont,1) = SIG(i);
                bb(cont,1) = i;
            end
        end
    end
    if ~exist('aa','var'); aa = NaN; bb = NaN; end
elseif band == -1
    cont = 0;
    for i = 2:length(SIGband)-1
        if SIGband(i) <= SIGband(i-1) && SIGband(i) < SIGband(i+1)
            if SIGband(i) > 0
                cont = cont+1;
                aa(cont,1) = SIG(i);
                bb(cont,1) = i;
            end
        end
    end
    if ~exist('aa','var'); aa = NaN; bb = 1; end
end
