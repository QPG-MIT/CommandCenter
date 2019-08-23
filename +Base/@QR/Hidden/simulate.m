cameraROI = [50,50];
cameraRes = [1000, 1000];
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
loc = [15, 9];
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
code = encode(30,24);
f = figure; ax = axes('parent',f); hold(ax,'on');
scatter(bits(:,1),bits(:,2)); hold on;
scatter(primary(:,1),primary(:,2));
scatter(secondary(:,1),secondary(:,2));

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