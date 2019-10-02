% Josh Popp

% 4.1 DISPLAY THE ANNULUS
% The location of each pixel
nX = 200; % number of columns
nY = 200; % number of rows
xj = 1 : nX; % x location of each column
yi = 1 : nY; % y location of each row
[xij,yij] = meshgrid(xj,yi); % x,y location of each pixel
% define an image of an annulus
cx = 75; % x component of center
cy = 75; % y component of center
r1 = 30; % inner radius in pixels
r2 = 45; % outer radius in pixels
% define the image using binary operations
I=((xij - cx).^2 + (yij - cy).^2 <= r2^2 )-((xij - cx).^2 + (yij - cy).^2 < r1^2 );

% display it
figure;
imagesc(I);
axis image;
title("I(x)");
grid on
set(gca,"ydir","normal"); % put origin at bottom left
hold on;

% 4.2 RECORD LANDMARKS
% Start at 3 o'clock and work counterclockwise
X = [111 68; 92 107; 56 107; 39 68; 56 41; 92 41];
Y = [180 120; 150 172; 90 172; 60 120; 90 63; 150 63];

% 4.3 FIND SCALE FACTOR 
s = sum(sum(Y)) / sum(sum(X));

% 4.4 TRANSFORM LANDMARKS
sX = X .* s;
figure;
scatter(X(:,1),X(:,2),'c', 'filled');
hold on;
scatter(sX(:,1),sX(:,2),'b', 'filled');
scatter(Y(:,1),Y(:,2),'r', 'filled');

% 4.5 NAIVE TRANSFORMATION
% initialize an image of all zeros
ITransformed = zeros(size(I));
for i = 2 : nY % loop through each row
    for j = 2 : nX % loop through each column
        % we are looking for the value to assign to Isx(j,i)
        % find the position to look at in the image J
        iLook = i*s;  % the inverse
        jLook = j*s;
        % round them to the nearest integer
        iLookRound = round(iLook);
        jLookRound = round(jLook);
        % check if we?re out of bounds,
        if iLookRound < 1 || iLookRound > nY || jLookRound < 1 || jLookRound > nX
            % if so, fill the image with the value zero
            ITransformed(j,i) = 0;
            else
            % otherwise, assign the value in our image at this point
            ITransformed(j,i) = I(jLookRound,iLookRound);
        end
        % don?t forget to index your images by (row,column) and not (x,y) !
    end
end

% display the new guy
figure;
imagesc(ITransformed);
axis image;
title("I'(x)");
grid on
set(gca,"ydir","normal"); % put origin at bottom left

% 4.6 OBSERER EQN TRANSFORMATION
% initialize an image of all zeros
ITransformed = zeros(size(I));
for i = 2 : nY % loop through each row
    for j = 2 : nX % loop through each column
        % we are looking for the value to assign to Isx(j,i)
        % find the position to look at in the image J
        iLook = i/s;  % the inverse
        jLook = j/s;
        % round them to the nearest integer
        iLookRound = round(iLook);
        jLookRound = round(jLook);
        % check if we?re out of bounds,
        if iLookRound < 1 || iLookRound > nY || jLookRound < 1 || jLookRound > nX
            % if so, fill the image with the value zero
            ITransformed(j,i) = 0;
            else
            % otherwise, assign the value in our image at this point
            ITransformed(j,i) = I(jLookRound,iLookRound);
        end
        % don't forget to index your images by (row,column) and not (x,y) !
    end
end

% display the new guy
figure;
imagesc(ITransformed);
axis image;
title("I'(x)");
grid on
set(gca,"ydir","normal"); % put origin at bottom left

hold on;
scatter(Y(:,1),Y(:,2),'r', 'filled');
