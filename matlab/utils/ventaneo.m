function [Nvent, M, iv, fv, wincleantot, wincleanEW, wincleanNS, wincleanVE, STALTAEW, STALTANS, STALTAVE] = ventaneo(porctrasl, ptosvent, Nn, EW, NS, VE, dt, tSTA, tLTA, Smax, Smin )
    Ntras = floor(porctrasl/100*ptosvent);
    iv = (1:ptosvent-Ntras:Nn).';
    fv = iv+ptosvent-1;
    % cambmax = [];
    cambmax = find(fv>Nn);
    if ~isempty(cambmax)
        fv(cambmax) = Nn; 
    end
    tvent = [iv,fv,fv-iv+ones(length(iv),1)];
    elim = find(and(fv==Nn,tvent(:,3)<ptosvent));
    % if length(elim) >= 2; tvent(elim(2:end),:) = []; iv(elim(2:end)) = []; fv(elim(2:end)) = []; end
    tvent(elim,:) = []; 
    iv(elim) = []; 
    fv(elim) = [];
    M = length(iv);

    % ELIMINA LAS VENTANAS MÁS ENERGÉTICAS DE LA SEÑAL EN SEGUNDOS
    if Smax == 0
        wincleantot = ones(M,1);
    else
        % wincleantot = zeros(M,1);
        [wincleanEW,STALTAEW] = picossig6(EW,dt,iv,fv,tSTA,tLTA,Smax,Smin);
        [wincleanNS,STALTANS] = picossig6(NS,dt,iv,fv,tSTA,tLTA,Smax,Smin);
        [wincleanVE,STALTAVE] = picossig6(VE,dt,iv,fv,tSTA,tLTA,Smax,Smin);
        wincleantot = wincleanEW.*wincleanNS.*wincleanVE;
    end

    Nvent = sum(sum(wincleantot));
