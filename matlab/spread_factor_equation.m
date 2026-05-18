function f = spread_factor_equation()
%SPREAD_FACTOR_EQUATION  Fit and plot the spread-factor calibration curve for WSP.
%
%   Water Sensitive Paper stains are larger than the droplets that caused them.
%   This function fits a log model  f(x) = a + b*log(x)  to the ratio
%   (stain diameter / droplet diameter) for six calibration points, then
%   plots the fitted curve over the range 1–1500 µm.
%
%   The fitted coefficients are used in DROPLET_STATISTICS to convert stain
%   diameters to real droplet diameters via the polynomial approximation:
%       droplet_diameter = 0.53549306 * stain_diameter
%                        - 0.000084839 * stain_diameter^2
%
%   Requires the MATLAB Curve Fitting Toolbox.
%
%   Returns
%   -------
%   f : cfit object from the Curve Fitting Toolbox

stain_diameter = 100:100:600;
drop_diameter  = [59, 109, 155, 200, 243, 285];
spread_factor  = stain_diameter ./ drop_diameter;

fit_type = fittype('a+b*log(x)');
f = fit(stain_diameter', spread_factor', fit_type);

stain_diam = 1:1500;
sf         = f.a + f.b * log(stain_diam);
figure, plot(stain_diam, sf)
grid on
xlabel('Stain Diameter', 'fontsize', 16)
ylabel('Spread Factor', 'fontsize', 16)
title('Spread Factor on Water Sensitive Paper', 'fontsize', 16)
equation_label = sprintf('f(x) = %4.2f + %4.2f * log(x)', f.a, f.b);
legend(equation_label, 'Location', 'southeast', 'fontsize', 16)
set(gca, 'FontSize', 16);

return
