function H_corrected = multiplicar_funcion_transferencia(H_original, freq, station)



H_transferencia = abs(H_transferencia)/max(abs(H_transferencia));

for k = 1:size(H_original,2)
    H_corrected(:,k) = H_original(:,k).*H_transferencia;
end
