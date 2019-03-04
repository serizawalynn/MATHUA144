clc;
alpha=[0,0,0]; %initial position
m=5;           %mass kg
v=[3,4,5];           %velocity of projectile m/s, vector form
r=.1 ;          %radius of projectile m
R=287.058;           %gas constant of the atmosphere. account for composition and humidity. J/KgK
T=293.5;           %Temperature of atmosphere K
P=101325;           %absolute atmospheric pressure Pa
omega=[0,0,7.2921159* 10^-5];       %rotational spin of planet  rads/s
w=[-pi/4,0,0];     %rotational spin of object rads/s
g=[0,0,-9.81] ;%gravity
lambda=pi/3 ;% angle of elevation
CD=.47 ;          %drag coeff for sphere
wind=[0,0,0]; %wind vector (velocity)
controltime= ((- v(3)) + sqrt( v(3)^2 - (2*g*alpha(3))))/-g(3);%flight time

%unit vectors:
unitomega=omega/( sqrt(omega(1)^2 + omega(2)^2 + omega(3)^2));
unitw=w/sqrt((w(1)^2 + w(2)^2 + w(3)^2));
unitv=v/sqrt((v(1)^2 + v(2)^2 + v(3)^2));

%netvalues
omegaprime=sqrt(omega(1)^2 + omega(2)^2 + omega(3)^2); %rotational spin of earth w/o vector form
wprime=sqrt( w(1)^2 + w(2)^2 + w(3)^2);
vprime=sqrt(v(1)^2 + v(2)^2 + v(3)^2);

%derivations
roh=P/(R*T);
CL=r*wprime/vprime;
A=pi* r^2;

%Drag force
FD=.5*roh*(vprime^2)*CD*A*(-unitv);

%magnus force
FM=(.5*CL*roh*A*(vprime^2))*(cross(unitw,unitv));

%Coriolis acceleration
aCor=[ (2*omegaprime*(sin(lambda))*v(2)),((-2*omegaprime*(sin(lambda))*v(1)) + (2*omegaprime*(cos(lambda))*v(3))),(-2*omegaprime*(cos(lambda))*v(2))];

%wind force
FW=(.5*roh*(sqrt(wind(1)^2 +wind(2)^2 +wind(3)^2))*A);

%netacceleration
acc=(g+aCor)+ ((FD+FM+FW)/m);

%time
t=(  -(v(3)+sqrt((v(3)^2 - 2*g*alpha(1))))/g);

%result vector
syms k
delta=[   alpha(1)+ (v(1)*k) + (.5*(acc(1)*(k^2))) , alpha(2) + (v(2)*k) + (.5*(acc(2)*(k^2))) , alpha(3) + (v(3)*k) + (.5*acc(3)*(k^2))];
control=[alpha(1)+(v(1)*k),alpha(2)+(v(2)*k),(alpha(3)+(v(3)*k)+(.5*g(3)*(k^2)))];
xpos=(delta(1));
ypos=(delta(2));
zpos=(delta(3));
controlxpos=control(1);
controlypos=control(2);
controlzpos=control(3);
fplot3(xpos,ypos,zpos,[0 t])
hold on
fplot3(controlxpos,controlypos,controlzpos,[0 t])
hold off
legend('Experiment','Control')
grid on
zlim([0 inf])
xlabel('West/East (-i,i) (m)')
ylabel('North/South (j,-j) (m)')
zlabel('Down/Up (-k,k) (m)')

