function [EWv, NSv, VEv, fechahmsvent] = division_ventanas_tiempo(ESTR, ptosvent, Nventefec, Ndias, wincleantot, iv, fv)

    EWv = (zeros(ptosvent,Nventefec));
    NSv = (zeros(ptosvent,Nventefec));
    VEv = (zeros(ptosvent,Nventefec));
    fechahmsvent = cell(Nventefec,1);
    cont = 0;
    for p = 1:Ndias
        ind = find(wincleantot{p}==1);
        for kk = 1:length(ind)
            q = ind(kk);
            cont = cont+1;
            EWv(:,cont) = ESTR.EWrot{p}(iv{p}(q):fv{p}(q));
            NSv(:,cont) = ESTR.NSrot{p}(iv{p}(q):fv{p}(q));
            VEv(:,cont) = ESTR.VE{p}(iv{p}(q):fv{p}(q));
            fechahmsvent(cont,1) = {[ESTR.vecfechahms{p},'_',num2str(q)]};
        end
    end


%% Codigo anterior <<inicia>>

%    Ndias = length(vecfechahms);
%    Nvini = Ndias*M; %Nvini = M;
%    EWv = (zeros(ptosvent,Nvini)); %single
%    NSv = (zeros(ptosvent,Nvini)); %single
%    VEv = (zeros(ptosvent,Nvini)); %single
%    fechahmsvent = cell(Nvini,1);
%    cont = 0;
%    for j = 1:Ndias
%        for kk = 1:M
%            cont = cont+1;
%            EWv(:,cont) = EW(iv(kk):fv(kk),j);
%            NSv(:,cont) = NS(iv(kk):fv(kk),j);
%            VEv(:,cont) = VE(iv(kk):fv(kk),j);
%            fechahmsvent(cont,1) = {[vecfechahms{j},'_',num2str(kk)]};
%        end
%    end
%    wincleantotlinea = reshape(wincleantot,[1,Nvini]);
%    ventok = find(wincleantotlinea~=0);

%% <<fin>> codigo anterior
