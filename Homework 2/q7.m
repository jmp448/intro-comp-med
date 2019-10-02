%%%% NUMBER 7

%% 7.1 GENERATE A GRID
% The location of each pixel
nX = 200; % number of columns
nY = 200; % number of rows
xj = 1 : nX; % x location of each column
yi = 1 : nY; % y location of each row
[xij,yij] = meshgrid(xj,yi); % x,y location of each pixel

%% 7.2 GENERATE THE LANDMARKS
% Set landmarks
% Note that landmarks are stored in COLUMNS, not rows
% Rows are x and y, respectively
X = [125, 150; 100, 50];
% R = [cos, sin; -sin, cos]*30
theta=30;
R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
Y = R * X;

%% 7.3 PLOT LANDMARKS AND GRID
down = 10; % downsampling is important so you can see things clearly
xijdown = xij(1:down:end,1:down:end);
yijdown = yij(1:down:end,1:down:end);
% This is a trick to plot a grid easily.
% We actually plot a 3D surface,
% but view it directly from above so it looks 2D
surf(xijdown,yijdown,ones(size(xijdown)),'facecolor','none','edgecolor','k');
grid off
hold on;
scatter(X(1,:),X(2,:),'o','MarkerFaceColor', 'c', 'MarkerEdgeColor', 'c');
scatter(Y(1,:),Y(2,:),'o','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
view(0,90);

%% 7.4 CALCULATE AN OPTIMAL 2x2 MATRIX TRANSFORMATION
A = (Y*X)*inv(X*X);
AX = A * X;

%% 7.5 TRANSFORMED LANDMARKS
Axij = zeros(200,200);
Ayij = zeros(200,200);
for i =1:200
    for j = 1:200
        Axij(i,j) = A(1,:)*[i;j];
        Ayij(i,j) = A(2,:)*[i;j];
    end
end
Axijdown = Axij(1:down:end,1:down:end);
Ayijdown = Ayij(1:down:end,1:down:end);
surf(Axijdown,Ayijdown,ones(size(Axijdown)),'facecolor','none','edgecolor','k');

grid off
hold on;
scatter(X(1,:),X(2,:),'o','MarkerFaceColor', 'c', 'MarkerEdgeColor', 'c');
scatter(Y(1,:),Y(2,:),'o','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
scatter(AX(1,:),AX(2,:),'o','MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
view(0,90);

%% 7.6 JACOBIAN
syms x y
phi_x = [cosd(theta)*x-sind(theta)*y];
phi_y = [sind(theta)*x+cosd(theta)*y];
jac11 = gradient(phi_x, [x]);
jac12 = gradient(phi_x, [y]);
jac21 = gradient(phi_y, [x]);
jac22 = gradient(phi_y, [y]);
jac = [jac11, jac12; jac21, jac22];
% just to check...
% jacobian([0.866*x-0.5*y; 0.5*x+0.866*y], [x, y]);

% Calculating the determinant at each point
% The determinant is independent of x and y so we just need to calculate it
% once
detJ = det(jac);
detJ = double(detJ);
Axijdown = Axij(1:down:end,1:down:end);
Ayijdown = Ayij(1:down:end,1:down:end);
hold off;
surf(Axijdown, Ayijdown, ones(size(Axijdown))*detJ);
colorbar

%% 7.7 GAUSSIAN KERNEL TRANSFORMATION
sigma=50;

k11 = exp(-norm(X(:,1)-X(:,1))^2/(2*sigma^2));
k12 = exp(-norm(X(:,1)-X(:,2))^2/(2*sigma^2));
k21 = exp(-norm(X(:,2)-X(:,1))^2/(2*sigma^2));
k22 = exp(-norm(X(:,2)-X(:,2))^2/(2*sigma^2));

Khat = [k11, k12; k21, k22];
V1 = [Y(1,1)-X(1,1);Y(1,2)-X(1,2)];
P1 = inv(Khat)*V1;

V2 = [Y(2,1)-X(2,1);Y(2,2)-X(2,2)];
P2 = inv(Khat)*V2;

%% 7.8 PLOT GAUSSIAN TRANSFORMED LANDMARKS AND GRID
% Landmarks
Xtransform = zeros(2, 2);
Xtransform(1,:) = X(1,:).' + V1;
Xtransform(2,:) = X(2,:).' + V2;

k11 = exp(-norm(X(:,1)-X(:,1))^2/(2*sigma^2));
k12 = exp(-norm(X(:,1)-X(:,2))^2/(2*sigma^2));
k21 = exp(-norm(X(:,2)-X(:,1))^2/(2*sigma^2));
k22 = exp(-norm(X(:,2)-X(:,2))^2/(2*sigma^2));

Khat = [k11, k12; k21, k22];
vX = Khat * [P1, P2];
Xtransform = X + vX.';

% Grid
% initialize to identity, we will add the displacment v
phix = xij;
phiy = yij;
for i = 1 : nY
    for j = 1 : nX
        % add the displacement for each p(k) in the sum
        for k = 1 : size(X,2) % number of landmarks
            %Kij the kernel evaluated at (j,i) - X(k)
            Kij = exp(-norm([j,i]-X(:,k))^2/(2*sigma^2));
            phix(j,i) = phix(j,i) + Kij*P1(k); % ... add the x component for p(k)
            phiy(j,i) = phiy(j,i) + Kij*P2(k); % ... add the y component for p(k)
        end
    end
end

Phixdown = phix(1:down:end,1:down:end);
Phiydown = phiy(1:down:end,1:down:end);

surf(Phixdown,Phiydown,ones(size(Phixdown)),'facecolor','none','edgecolor','k');
grid off
hold on;
scatter(X(1,:),X(2,:),'o','MarkerFaceColor', 'c', 'MarkerEdgeColor', 'c');
scatter(Y(1,:),Y(2,:),'o','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
scatter(Xtransform(1,:),Xtransform(2,:),'o','MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
view(0,90);

%% 7.9 Jacobian of the Gaussian
hold off
detJ = zeros(size(xij));
for i = 1:nX
    for j = 1:nY
        w = [i;j];
        dphixdx = 1 + exp(-norm(w-X(:,1))^2/(2*sigma^2))*P1(1)*(-1/sigma^2)*norm(w-X(:,1))+exp(-norm(w-X(:,2))^2/(2*sigma^2))*P1(2)*(-1/sigma^2)*norm(w-X(:,2));
        dphixdy = exp(-norm(w-X(:,1))^2/(2*sigma^2))*P1(1)*(-1/sigma^2)*norm(w-X(:,1))+exp(-norm(w-X(:,2))^2/(2*sigma^2))*P1(2)*(-1/sigma^2)*norm(w-X(:,2));
        dphiydx = exp(-norm(w-X(:,1))^2/(2*sigma^2))*P2(1)*(-1/sigma^2)*norm(w-X(:,1))+exp(-norm(w-X(:,2))^2/(2*sigma^2))*P2(1)*(-1/sigma^2)*norm(w-X(:,2));
        dphiydy = 1 + exp(-norm(w-X(:,1))^2/(2*sigma^2))*P2(1)*(-1/sigma^2)*norm(w-X(:,1))+exp(-norm(w-X(:,2))^2/(2*sigma^2))*P2(1)*(-1/sigma^2)*norm(w-X(:,2));
        detJ(i,j) = dphixdx*dphiydy - dphixdy*dphiydx;
    end
end
detJdown = detJ(1:down:end, 1:down:end);
surf(Phixdown,Phiydown,detJdown)
colorbar
view(0,-90);
%% 7.10 Is this a freakin joke its 2am
R = [cosd(45), -sind(45); sind(45), cosd(45)];
Y = R * X;
        