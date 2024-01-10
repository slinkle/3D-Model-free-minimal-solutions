% generate points
function world=generate_random_points(totalNum,outlierRate,featureNoise,featureNoiseBound,cubesize,K,R,t)
% cubePose=load('cube_voteMatrix.txt');
dismatchNum = floor(totalNum*outlierRate/100.0);
matchNum=totalNum-dismatchNum;
% generate match 2d-3d points with noise

points2d=[];
bearing2d=[];
points3d=[];
plot3(t(1),t(2),t(3),'gdiamond');hold on;
for i=1:matchNum
    while true
        noise_p = randn(2,1)*featureNoise; % standard deviation is yaw_eb
        if norm(noise_p) <= featureNoiseBound
            break;
        end
    end
    while true
        or_points3d = (rand(3,1)-0.5).*cubesize; % [-cubesize,cubesize]
        pC0=R*or_points3d+t;
        if pC0(3)<1
            continue
        end
        pC0=pC0/pC0(3);
        p2d=K*pC0+[noise_p;0];
        if p2d(1)<0 ||p2d(1)>K(1,3)*2 ||p2d(2)<0 ||p2d(2)>K(2,3)*2
            continue
        end
        points2d=[points2d,p2d]; % generate noise observation
        bearing2d = [bearing2d,K\p2d];
        points3d=[points3d,or_points3d];
        break;
    end
    plot3(or_points3d(1),or_points3d(2),or_points3d(3),'ro');hold on;
end

% generate dismatch 2d-3d points with random pose
dispoints3d = (rand(3,dismatchNum)-0.5).*cubesize*10;
for i=1:dismatchNum
    disangle = (-pi + (pi-(-pi)).*randn(1,3));
    disR = eul2rotm(disangle);
    dist = (-5 + (5-(-5)).*randn(3,1));
    dispC0=disR*dispoints3d(:,i)+dist;
    dispC0=dispC0/dispC0(3);
    points2d = [points2d,K*dispC0];
    bearing2d = [bearing2d,dispC0];
    points3d = [points3d,dispoints3d(:,i)];
    plot3(dispoints3d(1,i),dispoints3d(2,i),dispoints3d(3,i),'k.');hold on;
end
axis equal;hold off;drawnow
world.points2d=points2d;
world.bearing2d=bearing2d;
world.points3d=points3d;
end