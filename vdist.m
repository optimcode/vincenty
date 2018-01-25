function [dists, A12, A21] =  vdist(lat1, lon1, lat2, lon2)
%   Calculates distances and bearings between points.
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
%     lat1, lon1: the initial point coodinates (in degrees) can be a vector
%     lat2, lon2: the final point coodinates (in degrees) can be a vector
%     accuracy: accuracy for the vincenty convergence (optional)
%
%   Returns:
%     s, A12, A21, distance (km), initial bearing (deg), and back bearing (deg).

accuracy=1.0e-12;
if ~isvector(lat1) || ~isvector(lat2) || ~isvector(lon1)|| ~isvector(lon2)
    error('input should be a vector');
end

if (numel(lat1) ~= numel(lat2)) || (numel(lon1) ~= numel(lon2)) || (numel(lon1) ~= numel(lat1)) 
    error('inputs should have the same lengths');
end

lat1 = lat1(:);
lat2 = lat2(:);
lon1 = lon1(:);
lon2 = lon2(:);

%%

a = 6378137.0;        % semi-major axis (km), WGS84
f = 1./298.257223563; % flattening of the ellipsoid, WGS84
b = (1-f)*a;          % semi-minor axis

phi1 = deg2rad(lat1);
L1   = deg2rad(lon1);
phi2 = deg2rad(lat2);
L2   = deg2rad(lon2);

U1 = atan((1-f)*tan(phi1));
U2 = atan((1-f)*tan(phi2));
L = L2 - L1;

lmbda = L;
lastlmbda = 9.E40;

while abs(lmbda - lastlmbda) > accuracy
    
    lastlmbda = lmbda;
    
    sin_sigma = ((cos(U2).*sin(lmbda)).^2.0 + (cos(U1).*sin(U2) - sin(U1).*cos(U2).*cos(lmbda)).^2.0).^0.5;
    cos_sigma = sin(U1).*sin(U2) + cos(U1).*cos(U2).*cos(lmbda);
    sigma = atan2(sin_sigma, cos_sigma);
    
    sin_alpha = (cos(U1).*cos(U2).*sin(lmbda))./sin(sigma);
    cossq_alpha = 1 - sin_alpha.^2.0;
    
    cos2sigma_m = cos(sigma) - (2.*sin(U1).*(sin(U2)./cossq_alpha));
    
    C = (f/16.)*cossq_alpha.*(4. + f*(4. - 3.*cossq_alpha));
    
    lmbda = (L + ((1.0 - C).*f).*sin_alpha.*(sigma + C.*sin_sigma.*(cos2sigma_m + C.*cos_sigma.* (-1. + 2.*(cos2sigma_m.^2.0)))));
end

usq = cossq_alpha*(a^2.0 - b^2.0)/b^2.0;
A = 1 + (usq/16384.).*(4096. + usq.*(-768. + usq.*(320. - 175.*usq)));
B = (usq/1024.).*(256. + usq.*(-128. + usq.*(74. - 47.*usq)));
dsigma = (B.*sin(sigma).*(cos2sigma_m + 0.25*B.*(cos(sigma).*(-1. + 2.*cos2sigma_m.^2.0) - (1./6.)*B.*cos2sigma_m.*(-3. + 4.*sin(sigma).^2.0).* (-3. + 4.0*cos2sigma_m.^2.0))));

dists = b*A.*(sigma-dsigma);

A12 = atan2(cos(U2).*sin(lmbda), (cos(U1).*sin(U2) - sin(U1).*cos(U2).*cos(lmbda)))  ;
A21 = atan2(cos(U1).*sin(lmbda), (-sin(U1).*cos(U2) + cos(U1).*sin(U2).*cos(lmbda)));

if A21 < pi
    A21 = A21 + pi;
else
    A21 = A21 - pi;
end

A12 = (A12 + 2.*pi); % (2.*pi)
A21 = (A21 + 2.*pi); % (2.*pi)

A12 = mod(rad2deg(A12), 360);
A21 = mod(rad2deg(A21), 360);

end