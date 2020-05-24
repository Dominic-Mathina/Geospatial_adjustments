
a= 6378137;
f = 1/298.26;
b= a-(a*f);
e =2*f - f^2;
x=5421846.7945
y=3524556.035
z = 264941.186

%%Closed Formula%%

lam = atan(y/x);

p = sqrt(x^2 + y^2);
thita = atan((z*a)/(p*b));
e2 = (a^2-b^2)/b^2;

num = z+(e2*b*(sin(thita)^3));
den = p-(e*a*(cos(thita)^3));

phi = atan(num/den);


v=a/sqrt(1-(e * (sin(phi))^2));
yyyy = z/(sin(phi))
h1 = (z/sin(phi)) -((v*(1-e)))
h2= (sqrt((x^2 + y^2))/cos(phi))-v
%%%%Collate
phiFinal = (phi*180)/pi
lambda = (lam*180)/pi

