% glat = (1.421388889*pi)/180;
% glong = (38.30611111*pi)/180;
% h = 1565.3;

glat = (1.421388889*pi)/180;
glong = (38.30611111*pi)/180;
h = 1565.3;
a= 6378137;
f = 1/298.257223;
b= a-(a*f);
e =2*f - f^2;

v=a/sqrt(1-(e * (sin(glat))^2));
%%%%%%%%%%%%%%Problems %%%%%
%%Rectangular Coordinates conversion %%
x=(v+h)*cos(glat)*cos(glong)
y=(v+h)*cos(glat)*sin(glong)
z = ((v*(1-e))+h)*sin(glat)

%%%Geodetic Azimuth%%

% Gazimuth = Aazimuth-(n*tan(phi))

%%%%%%%%%%%%
% dAzimuth = (Along-glong)*sin(glat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%








