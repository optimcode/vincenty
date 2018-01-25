function [lat2, lon2, A21] = vreckon(lat, lon, dist_km, bearing)
%   Computes the coordinates from a point towards a bearing at given distances.
%
%   This uses Vincenty formula with an accuracy parameter used
%   to set the convergence of lambda. The default 1E-12 corresponds
%   to approximately 0.06 mm.
%
%   The formulas and nomenclature are from Vincenty, 1975:
%     https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
%   See also:
%     https://en.wikipedia.org/wiki/Vincenty's_formulae
%
%   Inputs:
%     lat,lon: the initial point coordinates (in degrees),
%     dist_km: the distance of the target point (in km)
%     bearing: the bearing angle (in degrees)
%     accuracy: accuracy for the vincenty convergence (optional)
%
%   Returns:
%     lat2, lon2, A21 point latitudes, longitudes and reverse bearings, all in degrees.


%%
accuracy=1.0E-12;
a = 6378137.0;        % semi-major axis (km), WGS84
f = 1./298.257223563; % flattening of the ellipsoid, WGS84
b = (1-f)*a;          % semi-minor axis

phi1 = deg2rad(lat);
L1   = deg2rad(lon);
alpha1 = deg2rad(bearing);
s = dist_km;

U1 = atan((1-f)*tan(phi1));

sigma1 = atan2(tan(U1), cos(alpha1));

sinalpha = cos(U1)*sin(alpha1);
cossq_alpha = (1. - sinalpha^2.0);
usq = cossq_alpha*(a^2.0-b^2.0)/b^2.0;

A = 1 + usq/16384. * (4096. + usq*(-768 + usq*(320.-175.*usq)));
B = usq/1024.*(256. + usq*(-128. + usq*(74.-47.*usq)));

sigma = s/(b*A);
lastsigma = 1.E40;

while max(abs(sigma - lastsigma)) > accuracy
    lastsigma = sigma;
    twosigmam = 2.*sigma1 + sigma;
    dsigma = (B*sin(sigma).*(cos(twosigmam) + 0.25*B *(cos(sigma).*(-1. + 2.*cos(twosigmam).^2.0) - (1./6.)*B*cos(twosigmam).* (-3. + 4. *sin(sigma).^2.0).* (-3. + 4.*cos(twosigmam).^2.0))));
    sigma = s/(b*A) + dsigma;
end

num = sin(U1).*cos(sigma) + cos(U1).*sin(sigma).*cos(alpha1);
den = (1.-f)*(sinalpha.^2.0 + (sin(U1).*sin(sigma) - cos(U1).*cos(sigma).*cos(alpha1)).^2.0).^0.5;

lat2 = atan2(num, den);

num = sin(sigma).*sin(alpha1);
den = cos(U1).*cos(sigma) - sin(U1).*sin(sigma).*cos(alpha1);
lmbda = atan2(num, den);

C = (f/16.)*cossq_alpha*(4. + f*(4. - 3.*cossq_alpha));

L = (lmbda - (1. - C)*f*sinalpha * (sigma + C*sin(sigma).*(cos(twosigmam) + C*cos(sigma).* (-1. + 2.*cos(twosigmam).^2.0))));
lon2 = L + L1;

num = sinalpha;
den = -sin(U1).*sin(sigma) + cos(U1).*cos(sigma).*cos(alpha1);
A21 = atan2(num, den);
A21 = (A21 + 3.*pi); % (2.*pi)

lat2 =  rad2deg(lat2);
lon2 = rad2deg(lon2);
A21 = mod(rad2deg(A21), 360);
end
