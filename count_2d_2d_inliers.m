function [inliersNum,pinliers]=count_2d_2d_inliers(esR,est,norm_xis,norm_xjs,essentialTh)
    E=skew(est)*esR;
    pcost=norm_xis.'*E*norm_xjs;
    pcost=abs(diag(pcost));
    pinliers=find(pcost<=essentialTh);
    inliersNum=length(pinliers);
end