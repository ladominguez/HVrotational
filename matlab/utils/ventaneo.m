function [Nventefec, M, iv, fv, winclean, STALTA] = ventaneo(porctrasl, ptosvent, ESTR, dt, tSTA, tLTA, Smax, Smin, Ndias)

    Ntras = floor(porctrasl/100*ptosvent);
    iv = {};
    fv = {};
    M = 0;
    for p = 1:Ndias
        Nn = length(ESTR.VE{p});
        iv{p} = (1:ptosvent-Ntras:Nn).';
        fv{p} = iv{p}+ptosvent-1;
        elim = find(fv{p} > Nn);
        iv{p}(elim) = [];
        fv{p}(elim) = [];
        % rev = [iv{p} fv{p} fv{p}-iv{p}+1];
        M = M+length(iv{p});
    end

    %% <<INICIA>> codigo anterior - se puede borrar 
    %Ntras = floor(porctrasl/100*ptosvent);
    %iv = (1:ptosvent-Ntras:Nn).';
    %fv = iv+ptosvent-1;
    %% cambmax = [];
    %cambmax = find(fv>Nn);
    %if ~isempty(cambmax)
    %    fv(cambmax) = Nn; 
    %end
    %tvent = [iv,fv,fv-iv+ones(length(iv),1)];
    %elim = find(and(fv==Nn,tvent(:,3)<ptosvent));
    %% if length(elim) >= 2; tvent(elim(2:end),:) = []; iv(elim(2:end)) = []; fv(elim(2:end)) = []; end
    %tvent(elim,:) = []; 
    %iv(elim) = []; 
    %fv(elim) = [];
    %M = length(iv);

    %% <<FIN>>  codigo anterior

    % ELIMINA LAS VENTANAS MÁS ENERGÉTICAS DE LA SEÑAL EN SEGUNDOS

    winclean = [];
    if Smax == 0
        for p = 1:Ndias
            winclean.tot{p} = ones(length(iv{p}),1);
        end
    else
        Nventefec = 0;
        for p = 1:Ndias
            [winclean.EW{p},STALTA.EW{p}] = picossig6(ESTR.EWrot{p},dt,iv{p},fv{p},tSTA,tLTA,Smax,Smin);
            [winclean.NS{p},STALTA.NS{p}] = picossig6(ESTR.NSrot{p},dt,iv{p},fv{p},tSTA,tLTA,Smax,Smin);
            [winclean.VE{p},STALTA.VE{p}] = picossig6(ESTR.VE{p},dt,iv{p},fv{p},tSTA,tLTA,Smax,Smin);
            winclean.tot{p} = winclean.EW{p}.*winclean.NS{p}.*winclean.VE{p};
            Nventefec = Nventefec+sum(winclean.tot{p});
        end
    end

    %if Smax == 0
   %     wincleantot = ones(M,1);
   % else
   %     % wincleantot = zeros(M,1);
   %     [wincleanEW,STALTAEW] = picossig6(EW,dt,iv,fv,tSTA,tLTA,Smax,Smin);
   %     [wincleanNS,STALTANS] = picossig6(NS,dt,iv,fv,tSTA,tLTA,Smax,Smin);
   %     [wincleanVE,STALTAVE] = picossig6(VE,dt,iv,fv,tSTA,tLTA,Smax,Smin);
   %     wincleantot = wincleanEW.*wincleanNS.*wincleanVE;
   % end

   % Nvent = sum(sum(wincleantot));
