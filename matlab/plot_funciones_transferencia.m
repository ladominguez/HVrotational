clear all
close all

% check https://www.passcal.nmt.edu/response_plotter
CMG40T = readtable('./instrument_responses/GuralpCMG40TResults.csv'); 
A0_40T = 5.714046123e8;

Tmax = 500; 
Tmin = 1/50;

fmin = 1/Tmax;
fmax = 1/Tmin;

N = 51; 

freq = logspace( log10(fmin), log10(fmax), N);

H_transferencia = obtener_funcion_transferencia(freq, 'AZUL');
loglog(freq, H_transferencia,'linewidth',2)
hold on

H_transferencia = obtener_funcion_transferencia(freq, 'TOME');
loglog(freq, H_transferencia,'r','linewidth',2)

H_transferencia = obtener_funcion_transferencia(freq, 'CRIS');
loglog(freq, H_transferencia,'g','linewidth',2)

H_transferencia = obtener_funcion_transferencia(freq, 'PURU');
loglog(freq, H_transferencia,'k','linewidth',2)

loglog(CMG40T.freq, CMG40T.mag/800, 'ko')

legend('Sillicon Audio', 'Reftek 151B', 'Guralp 40T', 'Trillium', 'location', 'best')
 
xlabel('Frequency [Hz]')
ylim([0.3, 1.2])



