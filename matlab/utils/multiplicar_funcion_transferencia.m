function H_corrected = multiplicar_funcion_transferencia(H_original, freq, estac, fecha)
	    s = 1j * freq * 2 * pi;
	    H_transferencia = s.^2./((s.^2 + 0.14807*s + 0.010966).*(s.^2 + 355.38*s + 63165)); % Reftek 131
	    H_transferencia = abs(H_transferencia)/max(abs(H_transferencia));
	    
        for i = 1:size(H_original,2)
            H_corrected(:,i) = H_original(:,i).*H_transferencia;
        end
