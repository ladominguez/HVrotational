function H_corrected = multiplicar_funcion_transferencia(H_original, freq, estac, fecha)

	% 

        % Tipos de sensores
        % 1. AZUL - Sillicon Audio
        % 2. TAHU - Sillicon Audio
        % 3. LAPA - Sillicon Audio
        % 4. BARA - Reftek
        % 5. NARA - Guralp 40T
        % 6. TOME - Reftek
        % 7. HUIR - Reftek
        % 8. CONE - Reftek
        % 9. PURU - Nanometrics
        % 10. CRIS - Guralp 40T

	s = 1j * freq * 2 * pi;

        if strcmp(station, 'AZUL') || strcmp(station, 'TAHU') || strcmp(station, 'LAPA')
            % TODO - Sillicon Audio
			% A0 = 9.52654e7
			% Zeros 2
			% Z1 = 0
			z2 = -4.71239e+03
			% Poles 4
			p1 = -2.69674e-02
			p2 = -2.08916E+03  -5.59148e3i
			p3 = -2.08916E+03  +5.59148e3i 
			p4 = -1.25664E+04
			H_transferencia = s.*(s-z2)/((s-p1).*(s-p2).*(s-p3).*(s-p4))	

        elseif strcmp(station, 'NARA') || strcmp(station, 'CRIS')
            % TODO - Guralp 40T
	    	% A0 = 2304000
	    	% K = 2*400
	    	% Zeros 2
	    	% Z1 = 0
	    	% Z2 = 0
	    	% Poles 5
	    	p1 = -11.78e-3 + 11.78e-3i
	    	p2 = -11.78e-3 - 11.78e-3i
	    	p3 = -160
	    	p4 = -80
	    	p5 = -180
	    	H_transferencia = s.^2./((s-p1).*(s-p2).*(s-p3).*(s-p4).*(s-p5))

        elseif strcmp(station, 'BARA') || strcmp(station, 'TOME') || strcmp('HUIR') || strcmp('CONE')
            % Reftek 151B
            % A0 = 63165
            % K  = 2000 V/m/s (differential output)
            % Zeros 2
            % Z1 = 0
            % Z2 = 0
            % Poles 4
            % p1 = -0.07405 + 0.07405i
            % p2 = -0.07405 - 0.07405i
            % p3 = -177.72 + 177.72i
            % p4 = -177.72 - 177.72i
            c1 = 0.14807;
            c2 = 0.010966;
            c3 = 355.38;
            c4 = 63165;
	    	H_transferencia = s.^2./((s.^2 + c1*s + c2).*(s.^2 + c3*s + c4)); 
		elseif strcmp(station,'PURU')
	    	% Trilium compact 120s poshole 2 generation TC120-PH2
	    	% A0 = 4.34493e17
	    	% Zeros 6
	    	% Z1=Z2=0
	    	z3 = -392
	    	z4 = -1960
	    	z5 = -1490 + 1740i
	    	z6 = -1490 - 1740i
	    	% Poles 11
	    	p1 = -0.03691 + 0.03702i
	    	p2 = -0.03691 - 0.03702i
	    	p3 = -343
	    	p4 = -370 + 467i
	    	p5 = -370 - 467i
	    	p6 = -836 + 1522i
	    	p7 = -836 - 1522i
	    	p8 = -4900 + 4700i
	    	p9 = -4900 - 4700i
	    	p10 = 6900
	    	p11 = -15000
	    	H_transferencia = s.^2.*(s-z3).*(s-z4).*(s-z5)./((s-p1).*(s-p2).*(s-p3).*(s-p4).*(s*p5).*(s-p6).*(s-p7).*(s-p8).*(s-p9).*(s-p10).*(s-p11))
        else
            disp(['Sensor no encontrado'])
        end

	H_transferencia = abs(H_transferencia)/max(abs(H_transferencia));

    for i = 1:size(H_original,2)
        H_corrected(:,i) = H_original(:,i).*H_transferencia;
    end