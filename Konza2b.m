LT=[98.14006920 2.97814320 98.1326053 2.9562094 98.1625891 2.9852127 98.13443795 2.97348473];
L=LT'

A=[98.14006605 2.978137923 1 0
    2.978137923 -98.14006605 0 1
    98.13260212	2.956203984 1 0
2.956203984 -98.13260212	0 1
98.16258587	2.985207603 1 0
2.985207603 -98.16258587	 0 1
98.13443508	2.973479102 1 0
2.973479102 -98.13443508 0 1];

N=A'*A;
Qxx=inv(N);
d=A'*L;
X=Qxx*d


% R matrix

R =[0.999999148210918   -0.009721496375278
    0.009721496375278 0.999999148210918]

   T = [0.028991410194061
   0.000787511467934]

PRM2 =[98.13160542
    2.96477070]
RI = inv(R);
Xunes = RI*PRM2+T

%prm2  9813160.542	296477.070

%prm4 9813593.054	298254.175

%kca3b 9813942.920	297489.790

%kc6b 9813036.700	295469.030
Xp1=[0.002675
-0.021724999
0.006425
0.015260001
-0.001600001
0.00514
0.034075001
-0.0363];
Yp1=[-0.008625
0.012225
0.00525
0.01362
0.023316667
-0.00128
-0.02985
-0.015725];






