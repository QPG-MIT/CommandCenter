cameraROI = [60,60];
cameraRes = [500, 500]; % approx row (y), col (x)
sample_size = [0 0] + 2000; % um

assert(Base.QR.length==25,'Only setup to test QRs with length 25')
modSize = Base.QR.module_size;
numMods = 5;
[Y,X] = meshgrid(linspace(modSize*numMods,modSize,numMods),...
    linspace(modSize,modSize*numMods,numMods));
bitsBase = [X(:), Y(:)]+Base.QR.d-modSize/2;
ncols = sample_size(1)/Base.QR.spacing_between;
nrows = sample_size(2)/Base.QR.spacing_between;

primary = NaN(0,2);
secondary = NaN(0,2);
bits = NaN(0,2);
% Pick random location and generate QRs
loc = [randi(ncols+1)-1, randi(nrows+1)-1];
for i = 0:1 % row
    for j = 0:1 % col
        col = loc(1) + i; % x
        row = loc(2) + j; % y
        offset = [col, row]*Base.QR.spacing_between;
        [c,~] = Base.QR.BasicBlock();
        c = c + offset;
        primary = [primary; c(1:3,:)];
        secondary = [secondary; c(4:end,:)];
        code = logical(encode(row, col));
        bits = [bits; bitsBase(code,:) + offset];
    end
end
% Generate image
f = figure('visible','off');
ax = axes('parent',f); hold(ax,'on');
axis(ax,'image');
QRs = hgtransform('parent',ax);
for i = 1:size(primary,1)
    circle(QRs,primary(i,:),Base.QR.r);
end
for i = 1:size(secondary,1)
    circle(QRs,secondary(i,:),Base.QR.secondary_r);
end
for i = 1:size(bits,1)
    d = Base.QR.module_size;
    px = bits(i,1) - d/2;
    py = bits(i,2) - d/2;
    rectangle(QRs,'Position',[px py d d],...
        'FaceColor','k','linestyle','none');
end
xlim(ax,[-1 1]*cameraROI(1)/2);
ylim(ax,[-1 1]*cameraROI(2)/2);
x = ax.XLim;
y = ax.YLim;
% Transform (arbitrary rotation and slight translation/scaling)
translation = (loc+0.5).*Base.QR.spacing_between + Base.QR.spacing/2;
scale = rand(1)+0.5;
rot = rand(1)*2*pi;
mt = makehgtform('translate',-[translation 0]); % Move to 0,0
ms = makehgtform('scale',scale);
mr = makehgtform('zrotate',rot);
QRs.Matrix = mr*ms*mt;

% Take image
set([f,ax],'units','points');
ax.Position(3:4) = cameraRes;
f.Position = [0,0,cameraRes+ax.Position(1:2)+5];
axis(ax,'off');
F = getframe(ax);
im = rgb2gray(F.cdata); % 0 - 255
imNoisy = imnoise(im,'gaussian',0.01);
delete(f);
%f = figure; colormap(f,'gray');
%ax = axes('parent',f);
%imagesc(ax,x,y,imNoisy);
%axis(ax,'image');
clear im;
im.image = flipud(imNoisy); % Fix image x/y issue
im.ROI = [x; y];

[pos,readInfo,fdebug] = Base.QR.reader(im,'sensitivity',1.2,'debug',true);
mscale = 1/sqrt(sum(readInfo.tform.T(1,1)^2 + readInfo.tform.T(2,1)^2));
mrot = asin(readInfo.tform.T(1,2));
mtranslation = readInfo.tform.T(3,1:2);

% taken from python git repo (updated for matlab; including indexing from 1)
function code = encode(n,m)
% encode(row,col)
VERSION = 3;
%|  Version  |   Row   |  Col    | Checksum(mod8) |
%|  4 bits   |  8 bits | 8 bits  |     3 bits     |
ENCODING.pad = [0,4] + 1;         % Pad(1) bits = 1, 5
ENCODING.version_bits = 4;
ENCODING.row_bits = 8;
ENCODING.col_bits = 8;
ENCODING.checksum_bits = 3;
ENCODING.n = 25;

version = fliplr(de2bi(VERSION,ENCODING.version_bits));
row = fliplr(de2bi(n,ENCODING.row_bits));
col = fliplr(de2bi(m,ENCODING.col_bits));
code = [version row col];
% Create checksum
checksum = fliplr(de2bi(mod(sum(code),2^ENCODING.checksum_bits),ENCODING.checksum_bits));
code = [code checksum];
% Put in pad
for pos = ENCODING.pad
    code = [code(1:pos-1) 1 code(pos:end)];
end
assert(length(code) == 25,'Code not correct length');
end

function h = circle(ax,c,r)
% Specifying any property again in varargin will override defaults here
d = r*2;
px = c(1)-r;
py = c(2)-r;
h = rectangle(ax,'Position',[px py d d],'Curvature',[1,1],...
    'FaceColor','k','linestyle','none');
end