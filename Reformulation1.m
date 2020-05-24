%%%%%%%%% Author, O.Jotham %%%%%%%%
%%%%%%%%% Date: 03.10.2017 %%%%%%%%
%%%% Close Range Photogrammetry %%%
format long;
%%% Data is first imported from file via a ui interface in excel format
 
%[fileName, pathName, filterIndex] = uigetfile({'*.xlsx'}, 'Select file(s)', 'Multiselect', 'on');
%Ini = uiimport('Initial_Conditions.xlsx');
%Ref = uiimport('Ref_Photo_coords.xlsx');

%Ini = Ini.data;
%Ref = Ref.data;

Ref = [-13.8100   38.5250   -77.4250   39.3100
       -12.7700  -16.2910   -73.5940  -15.4990
        63.7540   49.8200    3.3280   50.4920
        51.5410  -10.5730   -3.8230   -9.8890
       -12.8900    9.2330  -73.8600   10.0070
        59.2900   13.6320    2.0610   14.0650
        32.3510   31.8180  -28.9760   32.5210
        25.9870    7.0120  -28.8190    7.7040];

Ini = [2594.79   2022.36    695.81    2480.0   3033.0
       2607.84   2080.31    206.06    378.0    377.0
       3310.35   1986.42    816.31    559.0    469.0
       3263.1    2129.39    245.96     NaN       NaN
       2606.0    2080.0    439.0       NaN       NaN
       3320.0    2107.0    488.0       NaN       NaN
       3015.0    2008.0    645.0       NaN       NaN
       3013.0    2188.0    424.0       NaN       NaN];


%% focal lengths of photo1 and photo2 respectively %%
f1 = 165.89;
f2 = 165.77;


%% Refined photo coordinates of the object-point on photo1 %%
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
for xx = 1:1;
%%%% rotational elements for photo 1
r11a= cos(phi1)*cos(kap1);
r12a= cos(om1)*sin(kap1)+sin(om1)*sin(phi1)*cos(kap1);
r13a= sin(om1)*sin(kap1)-cos(om1)*sin(phi1)*cos(kap1);
r21a=-cos(phi1)*sin(kap1);
r22a= cos(om1)*cos(kap1)-sin(om1)*sin(phi1)*sin(kap1);
r23a= sin(om1)*cos(kap1)+cos(om1)*sin(phi1)*sin(kap1);
r31a= sin(phi1);
r32a=-sin(om1)*cos(phi1);
r33a= cos(om1)*cos(phi1);

%%%% rotational elements for photo 2
r11b= cos(phi2)*cos(kap2);
r12b= cos(om2)*sin(kap2)+sin(om2)*sin(phi2)*cos(kap2);
r13b= sin(om2)*sin(kap2)-cos(om2)*sin(phi2)*cos(kap2);
r21b=-cos(phi2)*sin(kap2);
r22b= cos(om2)*cos(kap2)-sin(om2)*sin(phi2)*sin(kap2);
r23b= sin(om2)*cos(kap2)+cos(om2)*sin(phi2)*sin(kap2);
r31b= sin(phi2);
r32b=-sin(om2)*cos(phi2);
r33b= cos(om2)*cos(phi2);

%%%%% The partial derivatives of the rotation elements with respect to the 
%%%%% three components of rotaton.
%%===================== For Photo 1 ===================== 
dr11dom1 = 0;
dr12dom1 = -sin(om1)*sin(kap1)+cos(om1)*sin(phi1)*cos(kap1);
dr13dom1 = cos(om1)*sin(kap1)+sin(om1)*sin(phi1)*cos(kap1);
dr21dom1 = 0;
dr22dom1 = -sin(om1)*cos(kap1)-cos(om1)*sin(phi1)*sin(kap1);
dr23dom1 = cos(om1)*cos(kap1) -sin(om1)*sin(phi1)*sin(kap1);
dr31dom1 = 0;
dr32dom1 = -cos(om1)*cos(phi1);
dr33dom1 = -sin(om1)*cos(phi1);
dr11dphi1 = -sin(phi1)*cos(kap1);
dr12dphi1 = sin(om1)*cos(phi1)*cos(kap1);
dr13dphi1 = -cos(om1)*cos(phi1)*cos(kap1);
dr21dphi1 = sin(phi1)*sin(kap1);
dr22dphi1 = -sin(om1)*cos(phi1)*sin(kap1);
dr23dphi1 = cos(om1)*cos(phi1)*sin(kap1);
dr31dphi1 = cos(phi1);
dr32dphi1 = sin(om1)*sin(phi1);
dr33dphi1 = -cos(om1)*sin(phi1);
dr11dkap1 = -cos(phi1)*sin(kap1);
dr12dkap1 = cos(om1)*cos(kap1)-sin(om1)*sin(phi1)*sin(kap1);
dr13dkap1 = sin(om1)*cos(kap1)+cos(om1)*sin(phi1)*sin(kap1);
dr21dkap1 = -cos(phi1)*cos(kap1);
dr22dkap1 = -cos(om1)*sin(kap1)-sin(om1)*sin(phi1)*cos(kap1);
dr23dkap1 = -sin(om1)*sin(kap1)+cos(om1)*sin(phi1)*cos(kap1);
dr31dkap1 = 0;
dr32dkap1 = 0;
dr33dkap1 = 0;

%%================ For Photo 2 ========================%%

dr11dom2 = 0;
dr12dom2 = -sin(om2)*sin(kap2)+cos(om2)*sin(phi2)*cos(kap2);
dr13dom2 = cos(om2)*sin(kap2)+sin(om2)*sin(phi2)*cos(kap2);
dr21dom2 = 0;
dr22dom2 = -sin(om2)*cos(kap2)-cos(om2)*sin(phi2)*sin(kap2);
dr23dom2 = cos(om2)*cos(kap2) -sin(om2)*sin(phi2)*sin(kap2);
dr31dom2 = 0;
dr32dom2 = -cos(om2)*cos(phi2);
dr33dom2 = -sin(om2)*cos(phi2);
dr11dphi2 = -sin(phi2)*cos(kap2);
dr12dphi2 = sin(om2)*cos(phi2)*cos(kap2);
dr13dphi2 = -cos(om2)*cos(phi2)*cos(kap2);
dr21dphi2 = sin(phi2)*sin(kap2);
dr22dphi2 = -sin(om2)*cos(phi2)*sin(kap2);
dr23dphi2 = cos(om2)*cos(phi2)*sin(kap2);
dr31dphi2 = cos(phi2);
dr32dphi2 = sin(om2)*sin(phi2);
dr33dphi2 = -cos(om2)*sin(phi2);
dr11dkap2 = -cos(phi2)*sin(kap2);
dr12dkap2 = cos(om2)*cos(kap2)-sin(om2)*sin(phi2)*sin(kap2);
dr13dkap2 = sin(om2)*cos(kap2)+cos(om2)*sin(phi2)*sin(kap2);
dr21dkap2 = -cos(phi2)*cos(kap2);
dr22dkap2 = -cos(om2)*sin(kap2)-sin(om2)*sin(phi2)*cos(kap2);
dr23dkap2 = -sin(om2)*sin(kap2)+cos(om2)*sin(phi2)*cos(kap2);
dr31dkap2 = 0;
dr32dkap2 = 0;
dr33dkap2 = 0;

J = [1;1;1;1;2;2;2;2;3;3;3;3;4;4;4;4];

%%================ The Jacobian A- Marix ========================
for i = 1:nn/2; t = 1:nn;
	k = (Pn(i)+1)/2;
	h = (Pn(t)+1)/2 +4;

	%m = i(:,1:13);
	%o = P(m);
	r = Pn(i);
	i = 2*r-1;

	
	Xbar1 = r11a*(Xi(k)-Xo1) + r12a*(Zi(k)-Zo1) - r13a*(Yi(k)-Yo1); 
	Ybar1 = r21a*(Xi(k)-Xo1) + r22a*(Zi(k)-Zo1) - r23a*(Yi(k)-Yo1);
	Zbar1 = r31a*(Xi(k)-Xo1) + r32a*(Zi(k)-Zo1) - r33a*(Yi(k)-Yo1);

	Xbar2 = r11b*(Xi(k)-Xo2) + r12b*(Zi(k)-Zo2) - r13b*(Yi(k)-Yo2);
 	Ybar2 = r21b*(Xi(k)-Xo2) + r22b*(Zi(k)-Zo2) - r23b*(Yi(k)-Yo2);
 	Zbar2 = r31b*(Xi(k)-Xo2) + r32b*(Zi(k)-Zo2) - r33b*(Yi(k)-Yo2);

	%%%% dXo1 %%%%
	dCRb(i,1)    = -f1*(-r11a*(Zbar1) + r31a*(Xbar1))/(Zbar1)^2;
	dCRb(i+1,1)  = -f1*(-r21a*(Zbar1) + r31a*(Ybar1))/(Zbar1)^2;
	dCRb(i+2,1)  = 0;
	dCRb(i+3,1)  = 0;

	%%%% dYo1 %%%%
	dCRb(i,2)    = -f1*(r13a*(Zbar1) - r33a*(Xbar1))/(Zbar1)^2;
	dCRb(i+1,2)  = -f1*(r23a*(Zbar1) - r33a*(Ybar1))/(Zbar1)^2;
	dCRb(i+2,2)  = 0;
	dCRb(i+3,2)  = 0;

	%%%% dZo1 %%%%
	dCRb(i,3)    = -f1*(-r12a*(Zbar1) + r32a*(Xbar1))/(Zbar1)^2;
	dCRb(i+1,3)  = -f1*(-r22a*(Zbar1) + r32a*(Ybar1))/(Zbar1)^2;
	dCRb(i+2,3)  = 0;
	dCRb(i+3,3)  = 0;

	%%%% dOm1 %%%%
	dCRb(i,4)    = -f1*((dr11dom1*(Xi(k)-Xo1) +dr12dom1*(Zi(k)-Zo1) - dr13dom1*(Yi(k)-Yo1))*(Zbar1) -(dr32dom1*(Zi(k)-Zo1) - dr33dom1*(Yi(k)-Yo1))*Xbar1)/(Zbar1)^2;
	dCRb(i+1,4)  = -f1*((dr21dom1*(Xi(k)-Xo1) +dr22dom1*(Zi(k)-Zo1) - dr23dom1*(Yi(k)-Yo1))*(Zbar1) -(dr32dom1*(Zi(k)-Zo1) - dr33dom1*(Yi(k)-Yo1))*Ybar1)/(Zbar1)^2;
	dCRb(i+2,4)  = 0;
	dCRb(i+3,4)  = 0;

	%%%% dPhi1 %%%%
	dCRb(i,5)    = -f1*((dr11dphi1*(Xi(k)-Xo1) +dr12dphi1*(Zi(k)-Zo1) - dr13dphi1*(Yi(k)-Yo1))*(Zbar1) - (dr31dphi1*(Xi(k)-Xo1) +dr32dphi1*(Zi(k)-Zo1) - dr33dphi1*(Yi(k)-Yo1))*Xbar1)/Zbar1^2;
	dCRb(i+1,5)  = -f1*((dr21dphi1*(Xi(k)-Xo1) +dr22dphi1*(Zi(k)-Zo1) - dr23dphi1*(Yi(k)-Yo1))*(Zbar1) - (dr31dphi1*(Xi(k)-Xo1) +dr32dphi1*(Zi(k)-Zo1) - dr33dphi1*(Yi(k)-Yo1))*Ybar1)/Zbar1^2;
	dCRb(i+2,5)  = 0;
	dCRb(i+3,5)  = 0;

	%%%% dKap1 %%%%
	dCRb(i,6)    = -f1*((dr11dkap1*(Xi(k)-Xo1) +dr12dkap1*(Zi(k)-Zo1) - dr13dkap1*(Yi(k)-Yo1))*Zbar1)/Zbar1^2;
	dCRb(i+1,6)  = -f1*((dr21dkap1*(Xi(k)-Xo1) +dr22dkap1*(Zi(k)-Zo1) - dr23dkap1*(Yi(k)-Yo1))*Zbar1)/Zbar1^2;
	dCRb(i+2,6)  = 0;
	dCRb(i+3,6)  = 0;

	%%%% dXo2 %%%%
	dCRb(i,7)    = 0;
	dCRb(i+1,7)  = 0;
	dCRb(i+2,7)  = -f2*(-r11b*(Zbar2) + r31b*(Xbar2))/(Zbar2)^2;
	dCRb(i+3,7)  = -f2*(-r21b*(Zbar2) + r31b*(Ybar2))/(Zbar2)^2;

	%%%% dYo2 %%%%
	dCRb(i,8)    = 0; 
	dCRb(i+1,8)  = 0; 
	dCRb(i+2,8)  = -f2*(r13b*(Zbar2) - r33b*(Xbar2))/(Zbar2)^2;
	dCRb(i+3,8)  = -f2*(r23b*(Zbar2) - r33b*(Ybar2))/(Zbar2)^2;

	%%%% dZo2 %%%%
	dCRb(i,9)    = 0;
	dCRb(i+1,9)  = 0;
	dCRb(i+2,9)  = -f2*(-r12b*(Zbar2) + r32b*(Xbar2))/(Zbar2)^2;
	dCRb(i+3,9)  = -f2*(-r22b*(Zbar2) + r32b*(Ybar2))/(Zbar2)^2;	

	%%%% dOm2 %%%%
	dCRb(i,10)    = 0;
	dCRb(i+1,10)  = 0;
	dCRb(i+2,10)  = -f2*((dr11dom2*(Xi(k)-Xo2) +dr12dom2*(Zi(k)-Zo2) - dr13dom2*(Yi(k)-Yo2))*(Zbar2) - (dr32dom2*(Zi(k)-Zo2) - dr33dom2*(Yi(k)-Yo2))*Xbar2)/Zbar2^2;
	dCRb(i+3,10)  = -f2*((dr21dom2*(Xi(k)-Xo2) +dr22dom2*(Zi(k)-Zo2) - dr23dom2*(Yi(k)-Yo2))*(Zbar2) - (dr32dom2*(Zi(k)-Zo2) - dr33dom2*(Yi(k)-Yo2))*Ybar2)/Zbar2^2;

	%%%% dPhi2 %%%%
	dCRb(i,11)    = 0;
	dCRb(i+1,11)  = 0;
	dCRb(i+2,11)  = -f2*((dr11dphi2*(Xi(k)-Xo2) +dr12dphi2*(Zi(k)-Zo2) - dr13dphi2*(Yi(k)-Yo2))*(Zbar2) - (dr31dphi2*(Xi(k)-Xo2) + dr32dphi2*(Zi(k)-Zo2) - dr33dphi2*(Yi(k)-Yo2))*Xbar2)/Zbar2^2;
	dCRb(i+3,11)  = -f2*((dr21dphi2*(Xi(k)-Xo2) +dr22dphi2*(Zi(k)-Zo2) - dr23dphi2*(Yi(k)-Yo2))*(Zbar2) - (dr31dphi2*(Xi(k)-Xo2) + dr32dphi2*(Zi(k)-Zo2) - dr33dphi2*(Yi(k)-Yo2))*Ybar2)/Zbar2^2;

	%%%% dKap2 %%%%
	dCRb(i,12)    = 0;
	dCRb(i+1,12)  = 0;
	dCRb(i+2,12)  = -f2*((dr11dkap2*(Xi(k)-Xo2) +dr12dkap2*(Zi(k)-Zo2) - dr13dkap2*(Yi(k)-Yo2))*Zbar2)/Zbar2^2;
	dCRb(i+3,12)  = -f2*((dr21dkap2*(Xi(k)-Xo2) +dr22dkap2*(Zi(k)-Zo2) - dr23dkap2*(Yi(k)-Yo2))*Zbar2)/Zbar2^2;
end
for i = 1:nn/2; 
	k = (Pn(i)+1)/2;

	%m = i(:,1:13);
	%o = P(m);
	r = Pn(i);
	i = 2*r-1;

	Xbar1 = r11a*(Xi(k)-Xo1) + r12a*(Zi(k)-Zo1) - r13a*(Yi(k)-Yo1); 
	Ybar1 = r21a*(Xi(k)-Xo1) + r22a*(Zi(k)-Zo1) - r23a*(Yi(k)-Yo1);
	Zbar1 = r31a*(Xi(k)-Xo1) + r32a*(Zi(k)-Zo1) - r33a*(Yi(k)-Yo1);

	Xbar2 = r11b*(Xi(k)-Xo2) + r12b*(Zi(k)-Zo2) - r13b*(Yi(k)-Yo2);
 	Ybar2 = r21b*(Xi(k)-Xo2) + r22b*(Zi(k)-Zo2) - r23b*(Yi(k)-Yo2);
 	Zbar2 = r31b*(Xi(k)-Xo2) + r32b*(Zi(k)-Zo2) - r33b*(Yi(k)-Yo2);

    Lo(i,1)   = (-f1*Xbar1)/Zbar1; 
    Lo(i+1,1) = (-f1*Ybar1)/Zbar1;
    Lo(i+2,1) = (-f2*Xbar2)/Zbar2;
    Lo(i+3,1) = (-f2*Ybar2)/Zbar2;
end
   L   = [Ref(1,:) Ref(2,:) Ref(3,:) Ref(4,:)];
   L   = L';
   A   = dCRb;
   N   = A'*A;
   dL  = L-Lo;
   AL  = A'*dL;
   Sol = N\AL;

   Xo1     =  Xo1    + Sol(1,1);
   Yo1     =  Yo1    + Sol(2,1);
   Zo1     =  Zo1    + Sol(3,1);
   om1     =  om1    + Sol(4,1);
   phi1    =  phi1   + Sol(5,1);
   kap1    =  kap1   + Sol(6,1);
   Xo2     =  Xo2    + Sol(7,1);
   Yo2     =  Yo2    + Sol(8,1);
   Zo2     =  Zo2    + Sol(9,1);
   om2     =  om2    + Sol(10,1);
   phi2    =  phi2   + Sol(11,1);
   kap2    =  kap2   + Sol(12,1);

end
Soln = [Xo1;Yo1;Zo1;om1;phi1;kap1;Xo2;Yo2;Zo2;om2;phi2;kap2];