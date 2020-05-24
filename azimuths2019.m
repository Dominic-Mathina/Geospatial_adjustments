
Along = 104.9976833;
glong = 104.9939556;
glat = 39.99348333;


%  Gazimuth = Aazimuth-(n*tan(phi))

%%%%%%%%%%%%
dAzimuth = (Along-glong)*sin((glat*pi)/180)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%