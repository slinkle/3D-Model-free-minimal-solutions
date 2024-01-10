function Tij=compute_absolute_2p_pose(res1,res2,Rjj2,tjj2)
    % estimate result
    Tij=[];
    tij_norm=res1(1:3,4);
    tij2_norm=res2(1:3,4);
%     tij_norm=[res1(3);0;res1(4)];
%     tij2_norm=[res2(3);0;res2(4)];
    Rij=res1(1:3,1:3);
    Rij2=res2(1:3,1:3);
    tji_norm=-Rij.'*tij_norm;
    tj2i_norm=-Rij2.'*tij2_norm;
    simk=5;
    simA=zeros(3,1);
    simC=tjj2;
    v1=tji_norm;
    v2=Rjj2*tj2i_norm;
    if dot(v1,v2)==1
        disp('colinear, no result');
        return;
    end
    simB=simA+simk*v1;
    simD=simC+simk*v2;
    vCA=simA-simC;
    vCD=simD-simC;
    vAB=simB-simA;
    dist=dot(vCA,cross(v1,v2));
    if abs(dist)>1e-5
        disp('not coplanar, no result')
        return;
    end
    num=norm(cross(vCA,vCD))/norm(cross(vCD,vAB));
    rhoAO=num*norm(vAB);
    simO=simA+rhoAO*v1;
    tji_est=simO-simA;
    tij_est=-Rij*tji_est;
    Tij=[Rij,tij_est;0 0 0 1];
end