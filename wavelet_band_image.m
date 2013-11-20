% [img] = wavelet_band_image(y,levels)
% converts a cell array of wavelet bands into an image;
%
% y - cell array of wavelet bands
% { scaling, h1,v1,d1,h2,v2,d2, ...., hn,vn,dn }
% levels - number of levels in decomposition
%
% img - image wavelet decomposition structured as 2d image
%
% orientations - 
% ---------
% |   | v |
% |-------|
% | h | d |
% ---------
%   h=1   v=2   d=3



function [img] = wavelet_band_image(y,levels)
    imgSize=size(y{1}).*2^levels;
    img=zeros(imgSize);
    k=1;
    putSection(y{k},[1 1],imgSize./2^(levels));
    k=k+1;

    OrientOffsets=[1 0;0 1;1 1];

    for iLevel=1:levels
        sectionSize=imgSize./(2^(levels-iLevel+1));
        for iOrient=1:3
            sectionOffset=OrientOffsets(iOrient,:).*sectionSize+[1 1];
            putSection(y{k},sectionOffset,sectionSize);
            k=k+1;
        end
    end

    function putSection(y,offset,sz)
        img(offset(1):offset(1)+sz(1)-1,offset(2):offset(2)+sz(2)-1)=y;
    end
end