function align

close all;

% make a figure
% it will have a 2x2 layout
pos = [100 100 512 512];
hFig = figure('position',pos,'toolbar','none');


gapLR = 0.125;
gapUD = 0.125;
width = (1-4*gapLR)/2;
height = (1-4*gapUD)/2;


% panels
% frame doesn't seem to actually draw anything
% if I want to put something into the panel, the units change
% when I put things in a frame, gap is now half the size
hPBL = uipanel('position',[0 0 0.5 0.5],'title','ZY View');
hPBR = uipanel('position',[0.5 0 0.5 0.5],'title','XY View');
hPTL = uipanel('position',[0 0.5 0.5 0.5],'title','Tools');
hPTR = uipanel('position',[0.5 0.5 0.5 0.5],'title','XZ View');
gapLR = gapLR*2;
gapUD = gapUD*2;
width = (1-2*gapLR);
height = (1-2*gapUD);



% xy on bottom right
hXY = axes('parent',hPBR,'position',[gapLR,gapUD,width,height]);
axis image
% axis square;
xlabel X
ylabel Y
zlabel Z
box on
h1(1) = patch('faces',[],'vertices',[],'facecolor','c','edgecolor','none','markerfacecolor','c','facelighting','gouraud','tag','atlas');
h2(1) = patch('faces',[],'vertices',[],'facecolor','r','edgecolor','none','markerfacecolor','r','facelighting','gouraud','tag','target');
light('position',[0 0 1])

% on top right we do XZ
hXZ = axes('parent',hPTR,'position',[gapLR,gapUD,width,height]);
axis image
% axis square;
xlabel X
ylabel Y
zlabel Z
view(0,0)
box on
h1(2) = patch('faces',[],'vertices',[],'facecolor','c','edgecolor','none','markerfacecolor','c','facelighting','gouraud','tag','atlas');
h2(2) = patch('faces',[],'vertices',[],'facecolor','r','edgecolor','none','markerfacecolor','r','facelighting','gouraud','tag','target');
light('position',[0 -1 0])

% on bottom left we do ZY
hZY = axes('parent',hPBL,'position',[gapLR,gapUD,width,height]);
axis image
% axis square;
xlabel X
ylabel Y
zlabel Z
view(-90,0)
camroll(-90)
box on
h1(3) = patch('faces',[],'vertices',[],'facecolor','c','edgecolor','none','markerfacecolor','c','facelighting','gouraud','tag','atlas');
h2(3) = patch('faces',[],'vertices',[],'facecolor','r','edgecolor','none','markerfacecolor','r','facelighting','gouraud','tag','target');
light('position',[-1 0 0])


% now place ui elements
% I'll put buttons on a grid
gridLR = 0.1;
gridUD = 0.1;
hLoadAtlas = uicontrol('parent',hPTL,'style','pushbutton','string','load atlas','units','normalized','position',[gridLR*1,1-gridUD*2,gridLR*3,gridUD],'callback',{@loaddata,'atlas'});
hLoadTarget = uicontrol('parent',hPTL,'style','pushbutton','string','load target','units','normalized','position',[gridLR*4,1-gridUD*2,gridLR*3,gridUD],'callback',{@loaddata,'target'});

% in the top left I want to display things
hMatrix = uicontrol('parent',hPTL,'style','text','units','normalized','position',[0.5*gridLR gridUD 1-gridLR 1-4*gridUD],'string','','backgroundcolor','w','tag','displaymatrix','horizontalalignment','center');



% also save

% now let's make buttons for translation
A = eye(4);

stepT = 1;

LPos = [gapLR 0 gapLR*0.5 gapUD*0.5];
RPos = [gapLR+width 0 gapLR*0.5 gapUD*0.5];
DPos = [0 gapUD gapLR*0.5 gapUD*0.5];
UPos = [0 gapUD+width gapLR*0.5 gapUD*0.5];

hXYL = uicontrol('parent',hPBR,'style','pushbutton','string','L','units','normalized','position',LPos,'userdata',A,'callback',{@apply,[1,0,0,-stepT;0,1,0,0;0,0,1,0;0,0,0,1]});
hXYR = uicontrol('parent',hPBR,'style','pushbutton','string','R','units','normalized','position',RPos,'userdata',A,'callback',{@apply,[1,0,0,stepT;0,1,0,0;0,0,1,0;0,0,0,1]});

hXYD = uicontrol('parent',hPBR,'style','pushbutton','string','D','units','normalized','position',DPos,'userdata',A,'callback',{@apply,[1,0,0,0;0,1,0,-stepT;0,0,1,0;0,0,0,1]});
hXYU = uicontrol('parent',hPBR,'style','pushbutton','string','U','units','normalized','position',UPos,'userdata',A,'callback',{@apply,[1,0,0,0;0,1,0,stepT;0,0,1,0;0,0,0,1]});


hXZL = uicontrol('parent',hPTR,'style','pushbutton','string','L','units','normalized','position',LPos,'userdata',A,'callback',{@apply,[1,0,0,-stepT;0,1,0,0;0,0,1,0;0,0,0,1]});
hXZR = uicontrol('parent',hPTR,'style','pushbutton','string','R','units','normalized','position',RPos,'userdata',A,'callback',{@apply,[1,0,0,stepT;0,1,0,0;0,0,1,0;0,0,0,1]});

hXZD = uicontrol('parent',hPTR,'style','pushbutton','string','D','units','normalized','position',DPos,'userdata',A,'callback',{@apply,[1,0,0,0;0,1,0,0;0,0,1,-stepT;0,0,0,1]});
hXZU = uicontrol('parent',hPTR,'style','pushbutton','string','U','units','normalized','position',UPos,'userdata',A,'callback',{@apply,[1,0,0,0;0,1,0,0;0,0,1,stepT;0,0,0,1]});
 


hZYL = uicontrol('parent',hPBL,'style','pushbutton','string','L','units','normalized','position',LPos,'userdata',A,'callback',{@apply,[1,0,0,0;0,1,0,0;0,0,1,-stepT;0,0,0,1]});
hZYR = uicontrol('parent',hPBL,'style','pushbutton','string','R','units','normalized','position',RPos,'userdata',A,'callback',{@apply,[1,0,0,0;0,1,0,0;0,0,1,stepT;0,0,0,1]});

hZYD = uicontrol('parent',hPBL,'style','pushbutton','string','D','units','normalized','position',DPos,'userdata',A,'callback',{@apply,[1,0,0,0;0,1,0,-stepT;0,0,1,0;0,0,0,1]});
hZYU = uicontrol('parent',hPBL,'style','pushbutton','string','U','units','normalized','position',UPos,'userdata',A,'callback',{@apply,[1,0,0,0;0,1,0,stepT;0,0,1,0;0,0,0,1]});


% now rotation
stepR = 5*pi/180; % radians
CCWPos = [0 0 gapLR*0.5 gapUD*0.5];
CWPos = [gapLR+width gapUD+height gapLR*0.5 gapUD*0.5];

hXYCCW = uicontrol('parent',hPBR,'style','pushbutton','string','CCW','units','normalized','position',CCWPos,'userdata',A,'callback',{@apply,[cos(stepR),-sin(stepR),0,0;sin(stepR),cos(stepR),0,0;0,0,1,0;0,0,0,1]});
hXYCW = uicontrol('parent',hPBR,'style','pushbutton','string','CW','units','normalized','position',CWPos,'userdata',A,'callback',{@apply,[cos(-stepR),-sin(-stepR),0,0;sin(-stepR),cos(-stepR),0,0;0,0,1,0;0,0,0,1]});

hXZCCW = uicontrol('parent',hPTR,'style','pushbutton','string','CCW','units','normalized','position',CCWPos,'userdata',A,'callback',{@apply,[cos(stepR),0,-sin(stepR),0;0,1,0,0;sin(stepR),0,cos(stepR),0;0,0,0,1]});
hXZCW = uicontrol('parent',hPTR,'style','pushbutton','string','CW','units','normalized','position',CWPos,'userdata',A,'callback',{@apply,[cos(-stepR),0,-sin(-stepR),0;0,1,0,0;sin(-stepR),0,cos(-stepR),0;0,0,0,1]});


hZYCCW = uicontrol('parent',hPBL,'style','pushbutton','string','CCW','units','normalized','position',CCWPos,'userdata',A,'callback',{@apply,[1,0,0,0;0,cos(-stepR),-sin(-stepR),0;0,sin(-stepR),cos(-stepR),0;0,0,0,1]});
hZYCW = uicontrol('parent',hPBL,'style','pushbutton','string','CW','units','normalized','position',CWPos,'userdata',A,'callback',{@apply,[1,0,0,0;0,cos(stepR),-sin(stepR),0;0,sin(stepR),cos(stepR),0;0,0,0,1]});



function loaddata(source,callbackdata,tag)
DialogTitle = ['Select ' upper(tag(1)) tag(2:end)];
FilterSpec = [{'*.byu';'*.lmk';'*.crv'},{'Surface (*.byu)';'Landmarks (*.lmk)';'Curves (*.crv)'}];
[FileName,PathName,FilterIndex] = uigetfile(FilterSpec,DialogTitle);
fid = fopen([PathName filesep FileName]);
if FilterIndex == 1 % surface
    line = fgetl(fid);
    data = sscanf(line,'%d %d %d');
    nV = data(2);
    nF = data(3);
    line = fgetl(fid);
    v = fscanf(fid,'%f %f %f\n',[3 nV])';
    f = abs(fscanf(fid,'%d %d %d\n',[3 nF])');
    set(findobj('tag',tag),'faces',f,'vertices',v)
elseif FilterIndex == 2 % landmarks
    % to do
    % if it is a landmark file
    % check the first line
    line = fgetl(fid);
    if ~strcmp(line,'Landmarks-1.0')
        error('Expecting a Landmarks-1.0 file')
    end
    line = fgetl(fid);
    nLandmarks = str2num(line);
    v = zeros(nLandmarks,3);
    f = zeros(nLandmarks,1);
    for i = 1 : nLandmarks
        line = fgetl(fid);
        name{i} = line;
        line = fgetl(fid);
        v(i,:) = sscanf(line,'%f %f %f 1 1');
        f(i) = i;
    end
    
    if strcmp(tag,'atlas')
        set(findobj('tag',tag),'faces',f,'vertices',v,'facecolor','none','edgecolor','none','markerfacecolor','c','markeredgecolor','none','markersize',5,'marker','o')
    else
        set(findobj('tag',tag),'faces',f,'vertices',v,'facecolor','none','edgecolor','none','markerfacecolor','r','markeredgecolor','none','markersize',5,'marker','o')
    end
    
    % hack for daniel on july 12, 2016
    set(findobj('tag',tag),'facevertexcdata',cellfun(@(x)str2num(x),name)','facecolor','interp','edgecolor','interp','markerfacecolor','flat')
    
    
elseif FilterIndex == 3 % txt curves.
%     keyboard
    line = fgetl(fid);
    nCurves = sscanf(line,'%d');    
    f = [];
    v = [];
    for i = 1 : nCurves
        line = fgetl(fid);
        nPoints = sscanf(line,'%f');
        p = fscanf(fid,'%f %f %f\n',[3 nPoints])';
        if nPoints > size(f,2)
            f = padarray(f,[0,nPoints-size(f,2)],'replicate','post');
            f_ = 1:nPoints;
        elseif nPoints < size(f,2)
            f_ = [1:nPoints,nPoints*ones(1,size(f,2)-nPoints)];
        end
        f = [f;f_+size(v,1)];
        v = [v;p];
    end 
    f = [f,f(:,1)];
    if strcmp(tag,'atlas')
        set(findobj('tag',tag),'faces',f,'vertices',v,'facecolor','none','edgecolor','c')
    else
        set(findobj('tag',tag),'faces',f,'vertices',v,'facecolor','none','edgecolor','r')
    end
end
fclose(fid);

function apply(source,callbackdata,B)
A = get(source,'userdata');

objects = findobj('tag','atlas');
v = get(objects(1),'vertices');

% we want to rotate about the object's center, not the origin
% v0 = max(v,[],1)*0.5 + min(v,[],1)*0.5;
v0 = mean(v,1);
B = [eye(3),v0';0 0 0 1]*B*[eye(3),-v0';0 0 0 1];


Ainv = inv(A);
v = bsxfun(@plus,v*Ainv(1:3,1:3)',Ainv(1:3,4)');



% do some display
A = B*A;
if sum(isnan(A))
    A = eye(4);
end

% these are all passed by reference, so no need to update
% but I'd really rather update
set(findobj('style','pushbutton'),'userdata',A);

set(findobj('tag','atlas'),'vertices',bsxfun(@plus,v*A(1:3,1:3)',A(1:3,4)'));

matrixText = [
    sprintf('%+07.2f %+07.2f %+07.2f %+07.2f', A(1,1), A(1,2), A(1,3), A(1,4))
    sprintf('%+07.2f %+07.2f %+07.2f %+07.2f', A(2,1), A(2,2), A(2,3), A(2,4))
    sprintf('%+07.2f %+07.2f %+07.2f %+07.2f', A(3,1), A(3,2), A(3,3), A(3,4))
    sprintf('%+07.2f %+07.2f %+07.2f %+07.2f', A(4,1), A(4,2), A(4,3), A(4,4))
    ];
set(findobj('tag','displaymatrix'),'string',matrixText)
