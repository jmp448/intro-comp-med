%% 1.5
clc 
% load template
I = double(imread('0001_CC_Con.png')>0);
% load target
IPrime = double(imread('0003_CC_Alz.png')>0);
% implement gradient descent alg
[ID,b] = grad_descent_chantelle(I,IPrime,0.001,72);
%% 1.5 plots
figure;
subplot(2,3,1);
imshow(I);
title('Template');
colorbar;
subplot(2,3,2);
imshow(IPrime);
title('Target');
colorbar;
subplot(2,3,3);
imshow(ID);
title('Translated template');
colorbar;
subplot(2,3,4);
imshow(I-IPrime);
title('template - target');
colorbar;
subplot(2,3,5);
imshow(ID-IPrime);
title('translated template - target');
colorbar;

%% 1.6 Image matching with splines and Jacobian calculations

% 1.7 
clc
alpha = 20;
sigma = 0.01;
epsilon = 0.015;
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
title('Translated template - target ');
colorbar;

%%
