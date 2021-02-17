% compute the sigmas of the bi- Gaussian approximation of the PSF
% using equations of Bo Zhang (SPIE 2006)


% lambda_em = emission wavelength
% lambda_ex = excitation wavelength
% NA = numerical aperture of the objective lens
% n = refractive index of the sample medium
% microscope = string that determines microscope type: 'widefield',
% 'confocal' or 'nipkow'

function [sigma_xy,sigma_z] = sigma_PSF_BoZhang(lambda_em, lambda_ex, NA, n, microscope)

if isempty(lambda_ex)
    lambda_ex = lambda_em;
end

switch microscope
    case 'widefield'    % Widefield Microscopy
        sigma_xy = 0.225 * lambda_em / NA ;
        sigma_z   = 0.78 * n * lambda_em / (NA*NA) ;

    case {'confocal', 'nipkow'}   % Laser Scanning Confocal Microscopy and Spinning Disc Confocal Microscopy
        sigma_xy = 0.225 / NA * lambda_ex * lambda_em / sqrt( lambda_ex^2 + lambda_em^2 ) ;
        sigma_z =  0.78 * n / (NA^2) *  lambda_ex * lambda_em / sqrt( lambda_ex^2 + lambda_em^2 ) ;

    otherwise
        error(['microscope = ',microscope,' is not a valid option in sigma_PSF_BoZhang !']);
end


% For comparison, below are the expressions used by G. Danuser in the
% Thomann et al. paper :
% 
% sigma_xy = 0.21 * lambda / NA;
% sigma_z = 0.66 * lambda * n / (NA*NA);