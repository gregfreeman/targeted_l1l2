% [y] = wavelet_band_tree(img,levels)
% get cell array of wavelet band for each level and orientation
%
% img - image wavelet decomposition structured as 2d image
% levels - number of levels in decomposition
%
% y - cell array of wavelet bands
% { scaling, h1,v1,d1,h2,v2,d2, ...., hn,vn,dn }
%
% orientations - 
% ---------
% |   | v |
% |-------|
% | h | d |
% ---------
%   h=1   v=2   d=3


function [y] = wavelet_band_tree(img,levels)
    nbands=3*levels+1;
    imgSize=size(img);
    y=cell(nbands,1);
    k=1;
    y{k}=getSection([1 1],imgSize./2^(levels));
    k=k+1;

    OrientOffsets=[1 0;0 1;1 1];

    for iLevel=1:levels
        sectionSize=imgSize./(2^(levels-iLevel+1));
        for iOrient=1:3
            sectionOffset=OrientOffsets(iOrient,:).*sectionSize+[1 1];
            y{k}=getSection(sectionOffset,sectionSize);
            k=k+1;
        end
    end

    function y=getSection(offset,sz)
        y=img(offset(1):offset(1)+sz(1)-1,offset(2):offset(2)+sz(2)-1);
    end
end