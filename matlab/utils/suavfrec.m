function SIGsuav = suavfrec(SIG,f,fs)
df = f(2)-f(1);
[Nf,Nv] = size(SIG);

porc1 = 2^(-1/(2*fs));
porc2 = 2^(+1/(2*fs));
fini = round((f*porc1-f(1))/df)+1;
ffin = round((f*porc2-f(1))/df)+1;

fini(fini<1) = 1;
fini(fini>Nf) = Nf;
ffin(ffin<1) = 1;
ffin(ffin>Nf) = Nf;

Ndat = ffin-fini+1;
SIGsuav = zeros(Nf,Nv);
for i = 1:Nf
    SIGsuav(i,:) = sqrt(sum(SIG(fini(i):ffin(i),:).^2,1)./Ndat(i));
end
