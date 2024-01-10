% input: 2d (xj1 and xj2) in view j and 2d (xi1 and xi2) in view i
% output: the theta in Rij and up-to-a-scale tij %view j to view i
function res=compute_E_2p(xjs,xis)
% equation 1 from one point correspondence
xj1=xjs(:,1);xj2=xjs(:,2);
xi1=xis(:,1);xi2=xis(:,2);
ui1=xi1(1);vi1=xi1(2);
uj1=xj1(1);vj1=xj1(2);
% vi1*sin(theta-phi)+vi1*uj1*cos(theta-phi)+vj1*sin(phi)-ui1*vj1*cos(phi) % a1*x1+b1*x2+c1*x3+d1*x4=0
% equation 2 from one point correspondence
ui2=xi2(1);vi2=xi2(2);
uj2=xj2(1);vj2=xj2(2);
% vi2*sin(theta-phi)+vi2*uj2*cos(theta-phi)+vj2*sin(phi)-ui2*vj2*cos(phi)

% ground truth x=[x1;x2;x3;x4];
% gtx=[sin(theta-phi);cos(theta-phi);sin(phi);cos(phi)]

A = [vi1,vi1*uj1,vj1,-ui1*vj1;...
     vi2,vi2*uj2,vj2,-ui2*vj2]; % Ax=0
z=null(A);
xbar1=z(:,1);
xbar2=z(:,2);
% gtlamda2=(gtx(2)*xbar1(1)-gtx(1)*xbar1(2))/(xbar1(1)*xbar2(2)-xbar2(1)*xbar1(2))
% gtlamda1=(gtx(1)-gtlamda2*xbar2(1))/xbar1(1)

% coefficients
a1 = xbar1(1)^2+xbar1(2)^2; b1=xbar2(1)^2+xbar2(2)^2; c1=2*xbar1(1)*xbar2(1)+2*xbar1(2)*xbar2(2); %equation 1
a3 = xbar1(3)^2+xbar1(4)^2; b3=xbar2(3)^2+xbar2(4)^2; c3=2*xbar1(3)*xbar2(3)+2*xbar1(4)*xbar2(4); %equation 2

% % test==1
% a1*gtlamda1^2+b1*gtlamda2^2+c1*gtlamda1*gtlamda2
% a3*gtlamda1^2+b3*gtlamda2^2+c3*gtlamda1*gtlamda2

% equation3=equation1-equation2
a2=a1-a3;b2=b1-b3;c2=c1-c3;
% %test==0
% a2*gtlamda1^2+b2*gtlamda2^2+c2*gtlamda1*gtlamda2 
a4=c2/(a1*c2-c1*a2);
b4=(c1*b2-c2*b1)/(a1*c2-c1*a2);
% %test==0
% gtlamda1^2-a4-b4*gtlamda2^2
a5=a1*b4+b1;
b5=a1*a4-1;

% (a5*gtlamda2^2+b5)^2
% c1^2*gtlamda2^2*(a4+b4*gtlamda2^2)

a6=a5^2-c1^2*b4;
b6=2*a5*b5-c1^2*a4;
c6=b5^2;
%test==0
% a6*gtlamda2^4+b6*gtlamda2^2+c6
res={};
delta=b6^2-4*a6*c6;
if delta<0
    return;
end
s_delta=sqrt(delta);

lamda2_squares=[];
la_s=(-b6+s_delta)/(2*a6);
if la_s>=0
    lamda2_squares=[lamda2_squares,la_s];
end
la_s=(-b6-s_delta)/(2*a6);
if la_s>=0
    [val,~]=min(abs(lamda2_squares-la_s));
    if val>0.0001
        lamda2_squares=[lamda2_squares,la_s];
    end
end
% lamda2_squares
lamdas=[];
for i=1:length(lamda2_squares)
    la_s=lamda2_squares(i);
    la2_abs=sqrt(la_s);
    la1_s=a4+b4*la_s;
    if la1_s>=0
        la1_abs=sqrt(la1_s);
        lamdas=[lamdas,[la1_abs;la2_abs],[la1_abs;-la2_abs],...
                       [-la1_abs;la2_abs],[-la1_abs;-la2_abs]];
    end
end
% lamdas
if isempty(lamdas)
    return;
end


%at most four solutions
for i=1:size(lamdas,2)
    lamda1=lamdas(1,i);
    lamda2=lamdas(2,i);
%     fprintf('lamda=%f,%f--------------\n',lamda1,lamda2)
    sin_theta_phi=lamda1*xbar1(1)+lamda2*xbar2(1);
    cos_theta_phi=lamda1*xbar1(2)+lamda2*xbar2(2);
    sin_phi=lamda1*xbar1(3)+lamda2*xbar2(3);
    cos_phi=lamda1*xbar1(4)+lamda2*xbar2(4);
    tr12=2*cos_theta_phi^3 - cos_theta_phi*(cos_phi^2 + sin_phi^2 + cos_theta_phi^2 + sin_theta_phi^2) + 2*cos_theta_phi*sin_theta_phi^2;
    tr21=cos_phi*(cos_phi^2 + sin_phi^2 + cos_theta_phi^2 + sin_theta_phi^2) - cos_phi*(2*cos_phi^2 + 2*sin_phi^2);
    tr23=sin_phi*(2*cos_phi^2 + 2*sin_phi^2) - sin_phi*(cos_phi^2 + sin_phi^2 + cos_theta_phi^2 + sin_theta_phi^2);
    tr32=-sin_theta_phi*(cos_phi^2 + sin_phi^2 + cos_theta_phi^2 + sin_theta_phi^2) + 2*sin_theta_phi^3 + 2*cos_theta_phi^2*sin_theta_phi;
    if sum(abs([tr12,tr21,tr23,tr32]))<0.00001
%         AA=[cos_phi,-sin_phi;sin_phi,cos_phi];
%         Ab=[sin_theta_phi;cos_theta_phi];
%         Atheta=AA\Ab;
%         sin_theta=Atheta(1);
%         cos_theta=Atheta(2);
        E=[0,cos_theta_phi,0;-cos_phi, 0, sin_phi; 0, sin_theta_phi, 0];%Eji
        exist_now=false;
        for j=1:size(res,2)
            resi=res{j};
            if norm(abs(resi)-abs(E))<0.0001
                exist_now=true;
            end
        end
        if ~exist_now
            res=[res,E];
        end
%         points1=xis(1:2,:).';
%         points2=xjs(1:2,:).';
%         [Rji, tji] = cv.recoverPose(E, points1, points2);
%         Rij=Rji.';
%         tij=-Rij*tji;
%         exist_or_not=false;
%         for j=1:length(res)
%             R=res{j}(1:3,1:3);
%             t=res{j}(1:3,4);
%             quatError=dcm2quat(R/Rij);%degree
%             rotaError=abs(2*atan(norm(quatError(1,1:3))/quatError(1,4)));
%             if rotaError>pi/2
%                 rotaError=pi-rotaError;
%             end
%             rotaError=rotaError/pi*180;
%             tranError=norm(t-tij);
%             if abs(rotaError)<0.001&&tranError<0.001
%                 exist_or_not=true;
%                 break;
%             end
%         end
%         if ~exist_or_not
%             res=[res,[Rij,tij]];
%         end
    end
end


