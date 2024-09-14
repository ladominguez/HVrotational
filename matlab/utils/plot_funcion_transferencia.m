function H = plot_funcion_transferencia(H_original):
	freq = logspace(-2, 2, 100) 
	K = 2000
	A0 = 63165
	s = 1j * freq * 2 * pi
	H_transferencia = s.^2./((s.^2 + 0.14807*s + 0.010966).*(s.^2 + 355.38*s + 63165)) % Reftek 131
	



