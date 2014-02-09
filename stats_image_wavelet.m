function y=stats_image_wavelet(img,levels,sz)

img=reshape(img,sz);
c=wavelet_band_tree(img,levels);
nbands=1+3*levels;
y(nbands,1)=struct();

for iBand=1:nbands
    x=c{iBand};
    n=numel(x);   
    
    [ y(iBand).alpha, y(iBand).sigma, y(iBand).mu] = GGDParameterEstimator(x(:));
    y(iBand).l2=norm(x(:),2);
    y(iBand).l1=norm(x(:),1);
    y(iBand).r=y(iBand).l1^2/y(iBand).l2^2/n;
    y(iBand).power=y(iBand).l2^2/n; 

end
