function H_corrected = multiplicar_funcion_transferencia(H_original, freq, station)



H_transferencia = obtener_funcion_transferencia(freq, station);

H_corrected = zeros(numel(H_transferencia), size(H_original,2));
for k = 1:size(H_original,2)
    H_corrected(:,k) = H_original(:,k).*abs(H_transferencia);
end
