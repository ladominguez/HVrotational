function [fLnorm,fVnorm,fTnorm,Lnorm,Vnorm,Tnorm,Nblanqfrec,facts] = ...
    F_normalizacionfrec(L,V,T,Nespec,band,onebit,dt,factap)

[~,Nv] = size(L);
NQ = Nespec/2+1;
f = linspace(0,1/(2*dt),NQ).';
Nblanqfrec = find(f>=0.1,1,'first');
facts = [0;0;0];

if onebit == 1
    L = sign(L);
    V = sign(V);
    T = sign(T);
end
fL = fft(L,Nespec);
fV = fft(V,Nespec);
fT = fft(T,Nespec);
fL = fL(1:NQ,:);
fV = fV(1:NQ,:);
fT = fT(1:NQ,:);

if band == 0
    fLnorm = fL;
    fVnorm = fV;
    fTnorm = fT;

elseif band == 2
    % División entre las suma de las energías potenciales de las 3 direcciones
    CEE = sqrt((abs(fL)).^2+(abs(fV)).^2+(abs(fT)).^2);
    CEE = fsuavi(abs(CEE),f,Nblanqfrec,0);
    fLnorm = fL./CEE;
    fVnorm = fV./CEE;
    fTnorm = fT./CEE;
    
elseif band == 3
    % Blanqueamiento espectral
    fLsmooth = fsuavi(abs(fL),f,Nblanqfrec,0);
    fVsmooth = fsuavi(abs(fV),f,Nblanqfrec,0);
    fTsmooth = fsuavi(abs(fT),f,Nblanqfrec,0);
    fLnorm = fL./fLsmooth;
    fVnorm = fV./fVsmooth;
    fTnorm = fT./fTsmooth;
end

% Transformada inversa de Fourier
if band == 0
    Lnorm = L;
    Vnorm = V;
    Tnorm = T;
elseif band == 2 || band == 3
    Lnorm = real(ifft([real(fLnorm(1,:));fLnorm(2:end-1,:); ...
        real(fLnorm(end,:));flipud(conj(fLnorm(2:end-1,:)))]));
    Vnorm = real(ifft([real(fVnorm(1,:));fVnorm(2:end-1,:); ...
        real(fVnorm(end,:));flipud(conj(fVnorm(2:end-1,:)))]));
    Tnorm = real(ifft([real(fTnorm(1,:));fTnorm(2:end-1,:); ...
        real(fTnorm(end,:));flipud(conj(fTnorm(2:end-1,:)))]));
    
    tap = repmat(tukeywin(Nespec,factap),1,Nv);
    Lnorm = Lnorm.*tap;
    Vnorm = Vnorm.*tap;
    Tnorm = Tnorm.*tap;
end
