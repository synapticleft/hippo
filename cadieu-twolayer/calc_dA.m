function dA = calc_dA(Z,I_E,m)

dA = bsxfun(@times,-m.I_noise_factors,I_E*Z.');%complex(I_E*real(Z).',I_E*imag(Z).'));%
