function [inliersNum,dpinliers]=count_3d_2d_inliers(esR,est,imagePoints,worldPoints,pointTh)
    pts = esR*worldPoints+repmat(est,1,size(worldPoints,2));
    testimagePoints = pts;
%     testimagePointsd=testimagePoints;
    testimagePoints = [testimagePoints(1,:)./testimagePoints(3,:);testimagePoints(2,:)./testimagePoints(3,:);ones(1,size(testimagePoints,2))];
    A=testimagePoints;
    vec = sum(A.^2);
    vec=sqrt(vec);
    [M,N] = size(A);
    B = repmat(vec,M,1);
    A = A./B;
    testimagePoints=A;
    pcost=testimagePoints.'*imagePoints;
    pcost=diag(pcost);
    pcost=1-pcost;
    %                 residue = imagePoints.'- testimagePoints.';
    %                 pcost = sum(residue.*residue,2);
    dpinliers=find(pcost<=pointTh);
    inliersNum=length(dpinliers);
end