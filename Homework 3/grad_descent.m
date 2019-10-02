function [ID,b] = grad_descent(I,IPrime,epsilon,nIter)
    % initialize b
    b = zeros(2,1);
    
for i = 1 : nIter
    % deform image. ID = I(x-b). output = 512x512
    ID = applyVToImage(I,b);
    
    % energy (cost function). C(b). output = a number
    E = 0.5*sum(sum(abs(ID-IPrime).^2));
    fprintf('Iteration %d of %d, energy is %g\n',i,nIter,E);  
    fprintf('b(1) is %d, b(2) is %d\n',b(1,1),b(2,1));
    
    % gradient of I. outputs are = 512x512
    [gradIx,gradIy] = gradient(ID);
    
    % gradient of cost function wrt b. 
    gradCostx = -sum(sum((ID-IPrime).*gradIx));
    gradCosty = -sum(sum((ID-IPrime).*gradIy));
    
    % update velocity. output = a number
    b(1,1) = b(1,1) - (gradCostx*epsilon);
    b(2,1) = b(2,1) - (gradCosty*epsilon);    
end

function ID = applyVToImage(I,b)
% deform image by composing with x-b (interpolating at specified points)
[X,Y] = meshgrid(1:size(I,2),1:size(I,1));
ID = interp2(I,X-b(1),Y-b(2),'linear',0);
