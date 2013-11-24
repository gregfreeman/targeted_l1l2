settings=[];

qmf = MakeONFilter('Symmlet',8);
sz=[512,512];
levels=5;
W_op=opWavelet(sz(1),sz(2),'Custom',qmf,levels);

for iImage=1:17
    settings.image=iImage;
    [image,settings]=CSEvaluationImageGetter(settings);
    x=W_op*image(:);
    z=1./normalize_image_wavelet(x,levels,sz);
    stats{iImage}=stats_image_wavelet(x.*z,levels,sz);
end

fields=fieldnames(stats{1}(1))

nImage=17;
for iBand=1:16    
    for iImage=1:nImage    
        for iField=1:length(fields)
            mu=0;
            for iImage2=1:nImage    
                if iImage~=iImage2            
                    mu=mu+stats{iImage2}(iBand).(fields{iField})/(nImage-1);                    
                end
            end
            stats_wo_image{iImage}(iBand).(fields{iField})=mu;
        end
    end
end



save l1l2_pristine_stats stats stats_wo_image
data.stats=stats;
data.stats_wo_image=stats_wo_image;
mat2json('l1l2_pristine_stats.json',data);