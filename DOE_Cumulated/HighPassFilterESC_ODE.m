function dEta = HighPassFilterESC_ODE(t, Eta, yin, wh)

% ODEs
dEta = -wh*Eta + wh*yin;

end