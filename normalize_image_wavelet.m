function y=normalize_image_wavelet(img,levels,sz)

img=reshape(img,sz);
c=wavelet_band_tree(img,levels);
d=cell(13,1);

for iBand=1:4
    d{iBand}=ones(size(c{iBand}));    
end
h=ones(5,5);
h=h./sum(h(:));
    
for iBand=5:length(c)
    x=c{iBand};
    % divisive normalization
    img_var=imfilter(x.^2,h,'same');
    img_std=sqrt(img_var);
    img_divisor=img_std+1;                
    d{iBand}=img_divisor;
end

y=wavelet_band_image(d,levels);
y=y(:);