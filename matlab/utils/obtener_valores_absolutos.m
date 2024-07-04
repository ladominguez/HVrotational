function [fNSvent, fEWvent, fVEvent, fHHvent]= obtener_valores_absolutos(fNSventnorm,fEWventnorm, fVEventnorm, ini, fin)

    fNSvent = abs(fNSventnorm(ini:fin,:));
    fEWvent = abs(fEWventnorm(ini:fin,:));
    fVEvent = abs(fVEventnorm(ini:fin,:));
    fHHvent = abs(sqrt((fNSventnorm(ini:fin,:).^2+fEWventnorm(ini:fin,:).^2)/2));
