clear all
close all

Tmax = 500; 
Tmin = 1/100;

fmin = 1/Tmax;
fmax = 1/Tmin;

N = 51; 

freq = logspace( log10(fmin), log10(fmax), N);

H_transferencia = obtener_funcion_transferencia(freq, 'AZUL');

loglog(freq, H_transferencia)


