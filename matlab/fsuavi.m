function S = fsuavi(SIG,f,Nsuav,fs)
% SIG = fHHventnorm;
% f = frec;
% Nsuav = 0;
% fs = 100;

[Nf,Nc] = size(SIG);
df = f(2)-f(1);
if Nsuav ~= 0 && ~rem(Nsuav,2); Nsuav = Nsuav+1; end
if Nsuav == 0
    fac1 = 2^-(1/(2*fs));
    fac2 = 2^(1/(2*fs));
    f1 = round(f*fac1/df)+1;
    f2 = round(f*fac2/df)+1;
    f1(f1<=1) = 1;
    f2(f2<1) = 1;
    f2(f2>Nf) = Nf;
    loc1 = find(f1==1);
    loc2 = loc1(end)+1;
    Nfprom = f2-f1+1;
    SIGsuma = cumsum(double(SIG).^2);
    S = (zeros(Nf,Nc)); %single
    S(loc1,:) = sqrt(SIGsuma(f2(loc1),:)./Nfprom(loc1));
    S(loc2:end,:) = sqrt((SIGsuma(f2(loc2:end,:),:)-SIGsuma(f1(loc2:end,:)-1,:))./Nfprom(loc2:end,:));
else
    pend = (1:2:Nsuav).';
    vec = [pend;Nsuav*ones(Nf-2*length(pend),1);flipud(pend)];
    f1 = (1:Nf).'-(vec-1)/2;
    f2 = (1:Nf).'+(vec-1)/2;
    loc1 = find(f1==1);
    loc2 = loc1(end)+1;
    Nfprom = f2-f1+1;
    SIGsuma = cumsum(double(SIG));
    S = (zeros(Nf,Nc)); %single
    S(loc1,:) = SIGsuma(f2(loc1),:)./Nfprom(loc1);
    S(loc2:end,:) = (SIGsuma(f2(loc2:end,:),:)-SIGsuma(f1(loc2:end,:)-1,:))./Nfprom(loc2:end,:);
end

% figure
% loglog(f,SIG(:,1)); hold on
% loglog(f,S(:,1)); hold on
% loglog(f,smooth(SIG(:,1),Nsuav),'--')
