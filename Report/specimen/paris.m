function Nf = paris(dsigma,a_init,a_end,C,m,beta)

% if beta is a function handle
try
    dKinvPowm = @(a) 1./(beta(a)*dsigma.*sqrt(pi.*a)).^m;
    Nf = 1/C*integral(dKinvPowm,a_init,a_end);
% if beta is constant
catch
    dKinvPowm = @(a) 1./(beta*dsigma.*sqrt(pi.*a)).^m;
    Nf = 1/C*integral(dKinvPowm,a_init,a_end);
end


end