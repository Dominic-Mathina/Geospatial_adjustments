%%%%%%%%% Author, O.Jotham %%%%%%%%
%%%%%%%%% Date: 03.10.2017 %%%%%%%%
%%%% Close Range Photogrammetry %%%

%%% Data is first imported from file via a ui interface in excel format
 
% // [fileName, pathName, filterIndex] = uigetfile({'*.xlsx'}, 'Select file(s)', 'Multiselect', 'on');
% // Ini = uiimport('Initial_Conditions.xlsx');
% // Ref = uiimport('Ref_Photo_coords.xlsx');

% // Ini = Ini.data;
% // Ref = Ref.data;

Ref = [-13.8100   38.5250  -77.4250   39.3100
       -12.7700  -16.2910  -73.5940  -15.4990
        63.7540   49.8200   3.3280   50.4920
        51.5410  -10.5730   -3.8230   -9.8890
       -12.8900    9.2330  -73.8600   10.0070
        59.2900   13.6320    2.0610   14.0650
        32.3510   31.8180  -28.9760   32.5210
        25.9870    7.0120  -28.8190    7.7040];

Ini = [2594.8    2022.4    695.8    2480.0   3033.0
       2607.8    2080.3    206.1    378.0    377.0
       3310.4    1986.4    816.3    559.0    469.0
       3263.1    2129.4    246.0       NaN       NaN
       2606.0    2080.0    439.0       NaN       NaN
       3320.0    2107.0    488.0       NaN       NaN
       3015.0    2008.0    645.0       NaN       NaN
       3013.0    2188.0    424.0       NaN       NaN];


%% focal lengths of photo1 and photo2 respectively %%
f1 = 165.89;
f2 = 165.77;
%% Rehined photo coordinates of the object-point on photo1 %%
Px1 = Ref(:,1);
Py1 = Ref(:,2);
%% Refined photo coordinates of the object-point on photo2 %%
Px2 = Ref(:,3);
Py2 = Ref(:,4);


%% Coordinates of the left exposure station %%
Xo1 = Ini(1,4);
Yo1 = Ini(3,4);
Zo1 = Ini(2,4);
%% Coordinates of the right exposure station %%
Xo2 = Ini(1,5);
Yo2 = Ini(3,5);
Zo2 = Ini(2,5);
%% Coordinates of control and new points %%
Xi   = Ini(:,1);
Yi   = Ini(:,2);
Zi   = Ini(:,3);


%% Initial values of rotational elements %%% 
om1  = 0;
phi1 = 0;
kap1 = 0;
om2  = 0;
phi2 = 0;
kap2 = 0;

n  = length(Px1);
nn = n*4;

p = 2*n-1;
q = 4;
Pn= repmat(1:2:p,[q 1]);
Pn=Pn(:);
P = 1:13;
P = P';
for xx = 1:3;
  %%%% rotational elements for photo 1
  r11a= cos(phi1)*cos(kap1);
  r12a= -cos(phi1)*sin(kap1);
  r13a= sin(phi1);
  r21a= cos(om1)*sin(kap1)+sin(om1)*sin(phi1)*cos(kap1);
  r22a= cos(om1)*cos(kap1)-sin(om1)*sin(phi1)*sin(kap1);
  r23a= -sin(om1)*cos(phi1);
  r31a= sin(om1)*sin(kap1)-cos(om1)*sin(phi1)*cos(kap1);
  r32a= sin(om1)*cos(kap1)+cos(om1)*sin(phi1)*sin(kap1);
  r33a= cos(om1)*cos(phi1);

  %%%% rotational elements for photo 2
  r11b= cos(phi2)*cos(kap2);
  r12b= -cos(phi2)*sin(kap2);
  r13b= sin(phi2);
  r21b= cos(om2)*sin(kap2)+sin(om2)*sin(phi2)*cos(kap2);
  r22b= cos(om2)*cos(kap2)-sin(om2)*sin(phi2)*sin(kap2);
  r23b= -sin(om2)*cos(phi2);
  r31b= sin(om2)*sin(kap2)-cos(om2)*sin(phi2)*cos(kap2);
  r32b= sin(om2)*cos(kap2)+cos(om2)*sin(phi2)*sin(kap2);
  r33b= cos(om2)*cos(phi2);



%%%%%%%%%%%%%%%%%%%%%% The Jacobian A-Matrix %%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nn; t = 1:nn/2;
	k = (Pn(i)+1)/2;
	h = (Pn(t)+1)/2 +4;

	%m = i(:,1:13);
	%o = P(m);
	r = Pn(i);
    
	i = 2*r-1;
	
	%j = 3*((Pn(t)+1)/2);
	%%%% dXo1 %%%%
	dCR(i,1)    = -f1*r11a - Px1(k)*r31a;
	dCR(i+1,1)  = -f1*r21a - Py1(k)*r31a;
	dCR(i+2,1)  = 0;
	dCR(i+3,1)  = 0;

	%%%% dYo1 %%%%
	dCR(i,2)    = f1*r13a + Px1(k)*r33a;
	dCR(i+1,2)  = f1*r23a + Py1(k)*r33a;
	dCR(i+2,2)  = 0;
	dCR(i+3,2)  = 0;

	%%%% dZo1 %%%%
	dCR(i,3)    = -f1*r12a - Px1(k)*r32a;
	dCR(i+1,3)  = -f1*r22a - Py1(k)*r32a;
	dCR(i+2,3)  = 0;
	dCR(i+3,3)  = 0;		

	%%%% dOm1 %%%%%
  dCR(i,4)    = f1*((-sin(om1)*sin(kap1)+cos(om1)*sin(phi1)*cos(kap1))*(Zi(k)-Zo1) - (cos(om1)*sin(kap1)+sin(om1)*sin(phi1)*cos(kap1))*(Yi(k)-Yo1)) -Px1(k)*((cos(om1)*sin(kap1)+sin(om1)*sin(phi1)*cos(kap1)) + (cos(om1)*cos(kap1) - sin(om1)*sin(phi1)*sin(kap1))*(Zi(k)-Zo1)+sin(om1)*cos(phi1)*(Yi(k)-Yo1));
  dCR(i+1,4)  = f1*((sin(om1)*sin(kap1)+cos(om1)*sin(phi1)*cos(kap1))*(Zi(k)-Zo1) + (-sin(om1)*cos(kap1)-cos(om1)*sin(phi1)*sin(kap1))*(Yi(k)-Yo1)) - Py1(k)*((cos(om1)*sin(kap1)+sin(om1)*sin(phi1)*cos(kap1)) + (cos(om1)*cos(kap1)-sin(om1)*sin(phi1)*sin(kap1))*(Zi(k)-Zo1)+sin(om1)*cos(phi1)*(Yi(k)-Yo1));
  dCR(i+2,4)  = 0;
  dCR(i+3,4)  = 0;

    %%%% dPhi1 %%%%%
  dCR(i,5)    = f1*(-sin(phi1)*cos(kap1)*(Xi(k)-Xo1) + sin(phi1)*sin(kap1)*(Zi(k)-Zo1) - cos(phi1)*(Yi(k)-Yo1)) + Px1(k)*(-cos(om1)*cos(phi1)*cos(kap1)*(Xi(k)-Xo1) + cos(om1)*cos(phi1)*sin(kap1)*(Zi(k)-Zo1) + cos(om1)*cos(phi1)*(Yi(k)-Yo1));
  dCR(i+1,5)  = f1*(sin(om1)*cos(phi1)*cos(kap1)*(Xi(k)-Xo1)- sin(om1)*cos(phi1)*sin(kap1)*(Zi(k)-Zo1)-sin(om1)*sin(phi1)*(Yi(k)-Yo1))+ Py1(k)*(cos(om1)*cos(phi1)*cos(kap1)*(Xi(k)-Xo1)+cos(om1)*cos(phi1)*sin(kap1)*(Zi(k)-Zo1)+cos(om1)*sin(phi1)*(Yi(k)-Yo1));
  dCR(i+2,5)  = 0;
  dCR(i+3,5)  = 0;

    %%%% dKap1 %%%%
  dCR(i,6)    = f1*(-cos(phi1)*sin(kap1)*(Xi(k)-Xo1) - cos(phi1)*cos(kap1)*(Zi(k)-Zo1)) + Px1(k)*((sin(om1)*cos(kap1)-sin(om1)*sin(phi1)*sin(kap1))*(Xi(k)-Xo1) + (-cos(om1)*sin(kap1)-sin(om1)*sin(phi1)*cos(kap1))*(Zi(k)-Zo1));
  dCR(i+1,6)  = f1*((cos(om1)*cos(kap1)-sin(om1)*sin(phi1)*sin(kap1))*(Xi(k)-Xo1)-(cos(om1)*sin(kap1)+sin(om1)*sin(phi1)*cos(kap1))*(Zi(k)-Zo1)) + Py1(k)*((sin(om1)*cos(kap1)-sin(om1)*sin(phi1)*sin(kap1))*(Xi(k)-Xo1) + (-cos(om1)*sin(kap1)-sin(om1)*sin(phi1)*cos(kap1))*(Zi(k)-Zo1));
  dCR(i+2,6)  = 0;
  dCR(i+3,6)  = 0;

	%%%% dXo2 %%%%
	dCR(i,7)    = 0;
	dCR(i+1,7)  = 0;
	dCR(i+2,7)  = -f2*r11b - Px2(k)*r31b;
	dCR(i+3,7)  = -f2*r21b - Px2(k)*r31b;

	%%%% dYo2 %%%%
	dCR(i,8)    = 0;
	dCR(i+1,8)  = 0;
	dCR(i+2,8)  = f2*r13b + Px2(k)*r33b;
	dCR(i+3,8)  = f2*r23b + Py2(k)*r33b;

	%%%% dZo2 %%%%
	dCR(i,9)    = 0;
	dCR(i+1,9)  = 0;
	dCR(i+2,9)  = -f2*r12b - Px2(k)*r32b;
	dCR(i+3,9)  = -f2*r22b - Py2(k)*r32b;		

    %%%% dOm2 %%%%
  dCR(i,10)   = 0;
  dCR(i+1,10) = 0;
  dCR(i+2,10) = f2*((-sin(om2)*sin(kap2)+cos(om2)*sin(phi2)*cos(kap2))*(Zi(k) - Zo2) - (cos(om2)*sin(kap2)+sin(om2)*sin(phi2)*cos(kap2))*(Yi(k)-Yo2)) -Px2(k)*((cos(om2)*sin(kap2)+sin(om2)*sin(phi2)*cos(kap2))*(Xi(k)-Xo2) + (cos(om2)*cos(kap2) - sin(om2)*sin(phi2)*sin(kap2))*(Zi(k)-Zo2)+sin(om2)*cos(phi2)*(Yi(k)-Yo1));
  dCR(i+3,10) = f2*((sin(om2)*sin(kap2)+cos(om2)*sin(phi2)*cos(kap2))*(Zi(k)-Zo2) + (-sin(om2)*cos(kap2)-cos(om2)*sin(phi2)*sin(kap2))*(Yi(k)-Yo2)) - Py2(k)*((cos(om2)*sin(kap2)+sin(om2)*sin(phi2)*cos(kap2)) + (cos(om2)*cos(kap2)-sin(om2)*sin(phi2)*sin(kap2))*(Zi(k)-Zo2)+sin(om2)*cos(phi2)*(Yi(k)-Yo2)); 
    
    %%%% dPhi2 %%%%
  dCR(i,11)   = 0;
  dCR(i+1,11) = 0;
  dCR(i+2,11) = f2*(-sin(phi2)*cos(kap2)*(Xi(k)-Xo2) + sin(phi2)*sin(kap2)*(Zi(k)-Zo2) - cos(phi2)*(Yi(k)-Yo2)) + Px2(k)*(-cos(om2)*cos(phi2)*cos(kap2)*(Xi(k)-Xo2) + cos(om2)*cos(phi2)*sin(kap2)*(Zi(k)-Zo2) + cos(om2)*cos(phi2)*(Yi(k)-Yo2));
  dCR(i+3,11) = f2*(sin(om2)*cos(phi2)*cos(kap2)*(Xi(k)-Xo2)- sin(om2)*cos(phi2)*sin(kap2)*(Zi(k)-Zo2)-sin(om2)*sin(phi2)*(Yi(k)-Yo2))+ Py2(k)*(cos(om2)*cos(phi2)*cos(kap2)*(Xi(k)-Xo2)+cos(om2)*cos(phi2)*sin(kap1)*(Zi(k)-Zo2)+cos(om1)*sin(phi1)*(Yi(k)-Yo2)); 

    %%%% dKap2 %%%%
  dCR(i,12)   = 0;
  dCR(i+1,12) = 0;
  dCR(i+2,12) = f2*(-cos(phi2)*sin(kap2)*(Xi(k)-Xo2) - cos(phi2)*cos(kap2)*(Zi(k)-Zo2)) + Px2(k)*((sin(om2)*cos(kap2)-sin(om2)*sin(phi2)*sin(kap2))*(Xi(k)-Xo2) + (-cos(om2)*sin(kap2)-sin(om2)*sin(phi2)*cos(kap2))*(Zi(k)-Zo2));
  dCR(i+3,12) = f2*((cos(om2)*cos(kap2)-sin(om2)*sin(phi2)*sin(kap2))*(Xi(k)-Xo2)-(cos(om2)*sin(kap2)+sin(om2)*sin(phi2)*cos(kap2))*(Zi(k)-Zo2)) + Py2(k)*((sin(om2)*cos(kap2)-sin(om1)*sin(phi2)*sin(kap2))*(Xi(k)-Xo2) + (-cos(om2)*sin(kap2)-sin(om2)*sin(phi2)*cos(kap2))*(Zi(k)-Zo2));
    
end


for i = 1:nn;
	k = (Pn(i)+1)/2;
	r = Pn(i);
	i = 2*r-1;
  
  Fo(i,1)   = f1*(r11a*(Xi(k)-Xo1) + r12a*(Zi(k)-Zo1) - r13a*(Yi(k)-Yo1)) + Px1(k)*(r31a*(Xi(k)-Xo1) + r32a*(Zi(k)-Zo1) - r33a*(Yi(k)-Yo1));
  Fo(i+1,1) = f1*(r21a*(Xi(k)-Xo1) + r22a*(Zi(k)-Zo1) - r23a*(Yi(k)-Yo1)) + Py1(k)*(r31a*(Xi(k)-Xo1) + r32a*(Zi(k)-Zo1) - r33a*(Yi(k)-Yo1));
  Fo(i+2,1) = f2*(r11b*(Xi(k)-Xo2) + r12b*(Zi(k)-Zo2) - r13b*(Yi(k)-Yo2)) + Px2(k)*(r31b*(Xi(k)-Xo2) + r32b*(Zi(k)-Zo2) - r33b*(Yi(k)-Yo2));
  Fo(i+3,1) = f2*(r21b*(Xi(k)-Xo2) + r22b*(Zi(k)-Zo2) - r23b*(Yi(k)-Yo2)) + Py2(k)*(r31b*(Xi(k)-Xo2) + r32b*(Zi(k)-Zo2) - r33b*(Yi(k)-Yo2));
end
  
  A  = dCR(1:16,1:12);
  N  = A'*A;
  L  = -Fo(1:16,1);
  dL = A'*L;
  Sol = N\dL;
  
  Xo1     = Xo1  + Sol(1,1);
  Yo1     = Yo1  + Sol(2,1);
  Zo1     = Zo1  + Sol(3,1);
  om1     = om1  + Sol(4,1);
  phi1    = phi1 + Sol(5,1);
  kap1    = kap1 + Sol(6,1);
  Xo2     = Xo2  + Sol(7,1);
  Yo2     = Yo2  + Sol(8,1);
  Zo2     = Zo2  + Sol(9,1);
  om2     = om2  + Sol(10,1);
  phi2    = phi2 + Sol(11,1);
  kap2    = kap2 + Sol(12,1);
  Xi(5,1) =  Sol(13,1);
  Yi(5,2) =  Sol(14,1);
  Zi(5,3) =  Sol(15,1);
  Xi(6,1) =  Sol(16,1);
  Yi(6,2) =  Sol(17,1);
  Zi(6,3) =  Sol(18,1);
  Xi(7,1) =  Sol(19,1);
  Yi(7,2) =  Sol(20,1);
  Zi(7,3) =  Sol(21,1);
  Xi(8,1) =  Sol(22,1);
  Yi(8,2) =  Sol(23,1);
  Zi(8,3) =  Sol(24,1);
end

Soln  = [Xo1;Yo1;Zo1;om1;phi1;kap1;Xo2;Yo2;Zo2;om2;phi2;kap2];