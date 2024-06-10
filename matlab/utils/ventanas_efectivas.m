function [EWv, NSv, VEv, fechahmsvent]=ventanas_efectivas(ESTR, ptosvent, Nventefec, Ndias, wincleantot, iv, fv)
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

%    EWv = EWv(:,ventok);
%    NSv = NSv(:,ventok);
%    VEv = VEv(:,ventok);
%    fechahmsvent = fechahmsvent(ventok);
    
%    ventNOefectiv = unique([find(sum(abs(EWv))==0) find(sum(abs(NSv))==0) find(sum(abs(VEv))==0)]);
%    ventefectiv = 1:Nvent;
%    ventefectiv(ventNOefectiv) = [];
%    Nv = length(ventefectiv);
%    EWv = EWv(:,ventefectiv);
%    NSv = NSv(:,ventefectiv);
%    VEv = VEv(:,ventefectiv);
%    fechahmsvent = fechahmsvent(ventefectiv);
