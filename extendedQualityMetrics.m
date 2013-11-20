
%% extendedQualityMetrics
function results=extendedQualityMetrics(results,inputImage,outputImage)
% results=extendedQualityMetrics(results,inputImage,outputImage)
%
% add extended (detailed bandwise metrics)
% the output results struct will contain fields for different metrics of
% image quality

results=testQualityMetrics(results,inputImage,outputImage);

x=results.W_op*outputImage(:);
x_true=results.W_op*inputImage(:);
x_int=results.W_op*results.intermediateImage(:);

levels=5;
sz=size(outputImage);
nbands=16;
z=1./normalize_image_wavelet(x,levels,sz);
z_true=1./normalize_image_wavelet(x_true,levels,sz);
z_int=1./normalize_image_wavelet(x_int,levels,sz);
results.norm_stat=stats_image_wavelet(x.*z,levels,sz);
results.norm_stat_true=stats_image_wavelet(x.*z_true,levels,sz);
results.norm_stat_int=stats_image_wavelet(x.*z_int,levels,sz);

results.normalized_alpha=zeros(nbands,1);
results.normalized_alpha_r=zeros(nbands,1);

x2d=reshape(x,sz);
imgpyr=wavelet_band_tree(x2d,levels);

for iBand=1:nbands
    imgband=imgpyr{iBand};
    h=ones(5,5);
    h=h./sum(h(:));
    img_mean=imfilter(imgband,h,'same');
    img_zero_mean=imgband-img_mean;
    img_var=imfilter(img_zero_mean.^2,h,'same');
    img_std=sqrt(img_var);
    img_std_max=max(img_std(:));
    results.img_std_max=img_std_max;
    img_divisor=img_std+1;
    img_normalized=img_zero_mean./img_divisor;
    % remove edge effects
    img_normalized2=img_normalized(4:end-3,4:end-3);
%     img_divisor2=img_divisor(4:end-3,4:end-3);
%     img_std2=img_std(4:end-3,4:end-3);
    % fit to generalized gaussian for testing normality
    [ results.normalized_alpha(iBand), ~, ~] = GGDParameterEstimator(img_normalized2(:));
    [ results.unnormalized_alpha, ~, ~] = GGDParameterEstimator(imgband(:));
    x_abs=abs(img_normalized2);
    x2=img_normalized2.^2;
    results.normalized_alpha_r(iBand)=sum(x_abs(:)/length(x2)).^2./(sum(x2(:)/length(x2))+1);

end