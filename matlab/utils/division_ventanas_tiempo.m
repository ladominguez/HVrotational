function [ventok, EWv, NSv, VEv, fechahmsvent] = division_ventanas_tiempo(vecfechahms, M, ptosvent, iv, wincleantot, fv, EW, NS, VE)

    Ndias = length(vecfechahms);
    Nvini = Ndias*M; %Nvini = M;
    EWv = (zeros(ptosvent,Nvini)); %single
    NSv = (zeros(ptosvent,Nvini)); %single
    VEv = (zeros(ptosvent,Nvini)); %single
    fechahmsvent = cell(Nvini,1);
    cont = 0;
    for j = 1:Ndias
        for kk = 1:M
            cont = cont+1;
            EWv(:,cont) = EW(iv(kk):fv(kk),j);
            NSv(:,cont) = NS(iv(kk):fv(kk),j);
            VEv(:,cont) = VE(iv(kk):fv(kk),j);
            fechahmsvent(cont,1) = {[vecfechahms{j},'_',num2str(kk)]};
        end
    end
    wincleantotlinea = reshape(wincleantot,[1,Nvini]);
    ventok = find(wincleantotlinea~=0);