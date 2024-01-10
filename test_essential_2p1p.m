
clear;
clc;

addpath('/home/slinkle/mexopencv');
% addpath('/home/yanmei/mexopencv');
% addpath('D:\mexopencv-2.4');
% set the camera intrinsic 904.9990358306457 904.9990358306457 686.1872787475586 498.9065246582031
K = [904.9990358306457, 0.000000e+00, 686.1872787475586;
     0.000000e+00, 904.9990358306457, 498.9065246582031;
     0.000000e+00, 0.000000e+00, 1.000000e+00];
 
theta = 30/180*pi;
phi = 12/180*pi;
rho=1.4;

theta2 = 8/180*pi;
phi2 = 45/180*pi;
rho2=2.5;

featureNoise=0;
% in view i (same as world)
% set a 3d point
% P=[0.1;0.1;0.5];
P=[0.06;0.13;0.5];
pi1=K*P;
pi1=pi1/pi1(3); % u v 1 the detected feature point
noise_p = randn(2,1)*featureNoise;
pi1=pi1+[noise_p;0];
xi1=inv(K)*pi1;% point on the normalized plane
% xi1=xi1./norm(xi1);

% in view j (transformed)
Rij=[cos(theta), 0, sin(theta);...
     0           1,      0    ;...
     -sin(theta),0, cos(theta)]; %Ry(theta) from j to i (from camera to world)
tij=rho*[sin(phi);0;cos(phi)];
sj=norm(tij);
tij_norm=tij/sj;
Rji=Rij.'; % E
tji=-Rji*tij; % E
tji_norm=tji./norm(tji);
Tij_norm=[Rij,tij_norm;0 0 0 1];
R=Rji; % E=skew(t)*R
t=tji;
% feature point
Pj1=Rji*P+tji;
pj1=K*(Rji*P+tji);
pj1=pj1/pj1(3);
noise_p = randn(2,1)*featureNoise;
pj1=pj1+[noise_p;0];
xj1=inv(K)*pj1;
% xj1=xj1./norm(xj1);

% test the epipolar constraint
Eji=skew(t)*R;
xj1.'*Eji*xi1
% E_plane = rho*[0,cos(theta-phi),0;-cos(phi),0,sin(phi);0,sin(theta-phi),0]
E=skew(tij)*Rij;
xi1.'*E*xj1

% theta2=theta+25/180*pi;
% Rij=[cos(theta2), 0, sin(theta2);...
%      0           1,      0    ;...
%      -sin(theta2),0, cos(theta2)]; %Ry(theta) from j to i (from camera to world)
% % rho=1.7;
% phi = 15/180*pi+30/180*pi;
% tij=rho*[sin(phi);0;cos(phi)];
% E=skew(tij)*Rij;
% xi1.'*E*xj1



% generate another point
P2=[0.2;0.08;0.4];
% in view i (same as world)
pi2=K*P2;
pi2=pi2/pi2(3); % u v 1 the detected feature point
noise_p = randn(2,1)*featureNoise;
pi2=pi2+[noise_p;0];
xi2=inv(K)*pi2;% point on the normalized plane
% in view j
Pj2=Rji*P2+tji;
pj2=K*(Rji*P2+tji);
pj2=pj2/pj2(3);
noise_p = randn(2,1)*featureNoise;
pj2=pj2+[noise_p;0];
xj2=inv(K)*pj2;
xi2.'*E*xj2

gt_res=[sin(theta);cos(theta);sin(phi);cos(phi)]
res=compute_Rt_2p([xj1,xj2],[xi1,xi2])
% we can get Rji and tji_norm


% IN VIEW J2
Rij2=[cos(theta2), 0, sin(theta2);...
     0           1,      0    ;...
     -sin(theta2),0, cos(theta2)]; %Ry(theta) from j to i (from camera to world)

tij2=rho2*[sin(phi2);0;cos(phi2)];
Rj2i=Rij2.'; % E
tj2i=-Rj2i*tij2; % E
Tij=[Rij,tij;0 0 0 1];
Tij2=[Rij2,tij2;0 0 0 1];
Tjj2=inv(Tij)*Tij2;
Rjj2=Tjj2(1:3,1:3);
tjj2=Tjj2(1:3,4);
Tj2j=inv(Tjj2);
Rj2j=Tj2j(1:3,1:3);
tj2j=Tj2j(1:3,4);

P=[0.26;0.17;0.3];
pi1=K*P;
pi1=pi1/pi1(3); % u v 1 the detected feature point
noise_p = randn(2,1)*featureNoise;
pi1=pi1+[noise_p;0];
xi1=inv(K)*pi1;% point on the normalized plane
% xi1=xi1./norm(xi1);
% feature point
Pj1=Rj2i*P+tj2i;
pj1=K*(Rj2i*P+tj2i);
pj1=pj1/pj1(3);
noise_p = randn(2,1)*featureNoise;
pj1=pj1+[noise_p;0];
xj1=inv(K)*pj1;

% compute sj using Rji tji_norm and xiq xj1 Tjj
R2=Rj2j*Rji;
t2A=Rj2j*tji_norm;
t2b=tj2j;
r1=mean([R2(1,1),R2(3,3)]);
r2=mean([R2(1,3),-R2(3,1)]);
x1=xi1(1);y1=xi1(2);
x2=xj1(1);y2=xj1(2);
A1=r2*x1*y2+y1-r1*y2;
A3=r1*x1*y2-x2*y1+r2*y2;
a1=t2A(1);
a3=t2A(3);
b1=t2b(1);
b3=t2b(3);
s=-(A1*b1+A3*b3)/(A1*a1+A3*a3)
sj
s1=tj2i(1);
s3=tj2i(3);
A1*s1+A3*s3



