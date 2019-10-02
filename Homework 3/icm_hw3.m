%% 1.5
clc 
% load template
I = double(imread('0001_CC_Con.png')>0);
% load target
IPrime = double(imread('0003_CC_Alz.png')>0);
% implement gradient descent alg
step = 0.02;
nIter = 200;

[ID,b] = grad_descent(I,IPrime,step,nIter);
%% 1.5 plots - translation results
figure;
subplot(2,2,1);
imshow(I);
title('Template');
colorbar;
subplot(2,2,2);
imshow(IPrime);
title('Target');
colorbar;
subplot(2,2,3);
imshow(I);
title('Template');
colorbar;
subplot(2,2,4);
imshow(I);
title('Translated Template');
colorbar;

%% 1.5 Plots - translation differences
figure;
subplot(1,2,1);
imshow(I - IPrime);
title('Template - Target');
colorbar;
subplot(1,2,2);
imshow(ID - IPrime);
title('Translation - Target');
colorbar;

%% 1.7 using splines for transformation
clc
alpha = 20;
sigma = 0.01;
epsilon = 0.0175;
nIter = 200;
% energy min = 604257 at e = 0.015, niter = 200;
% implementing splineimage()
[ID2,vx,vy] = splineImage(ID,IPrime,alpha,sigma,epsilon,nIter);

% 1.7 plots
figure;
subplot(1,3,1);
imshow(ID2);
colorbar;
title('Deformed image');
subplot(1,3,2);
imshow(IPrime);
title('Target');
colorbar;
subplot(1,3,3);
imshow(ID2-IPrime);
title('Deformed image - Target ');
colorbar;

%% 1.8 Jacobian
figure;
[dvxdx, dvxdy] = gradient(vx);
[dvydx, dvydy] = gradient(vy);
detj = zeros(size(vx,1), size(vx,2));
for i = 1:size(vx,1)
    for j = 1:size(vx,2)
        detj(i,j) = (1+dvxdx(i,j))*(1+dvydy(i,j))-(1+dvxdy(i,j))*(1+dvydx(i,j));
    end
end
down=10;
xi = 1 : size(vx,1); % x location of each column
yj = 1 : size(vx,2); % y location of each row
[xij, yij] = meshgrid(xi,yj);
xijdown = xij(1:down:end,1:down:end);
yijdown = yij(1:down:end,1:down:end);
detjdown = detj(1:down:end,1:down:end);
grid off
surf(xijdown,yijdown,detjdown);
title('Determinant of the Jacobian');
colorbar;

%%
figure;
down=10;
Idown = ID(1:down:end,1:down:end);
surf(xijdown,yijdown,Idown);
view(0,90);
