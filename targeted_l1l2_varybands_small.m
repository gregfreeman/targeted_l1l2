function targeted_l1l2_varybands

param1.field='delta';
% param1.values=num2cell(0.1:0.1:0.4);
param1.values=num2cell(0.05:0.05:0.95);

param2.field='image';
param2.values=num2cell(1:17);

param3.field='cheat_r';
param3.values={false};

param4.field='cheat_z';
param4.values={'em'};

param5.field='startBand';
param5.values={2,5,8,11,13};

param6.field='image_options';
param6.values={{'cropsize',[128 128]}};

paramset=[param1,param2,param3,param4,param5,param6];

events=struct();
events.runExperiment=@getExperiment;
events.loadInputData=@(settings) CSEvaluationImageGetter(settings);
events.evaluateMetrics= @extendedQualityMetrics;
events.storeOutputData= @pngStoreOutputData;
events.startTasks = @( foldername, ntasks ) startTasksHtmlPost( foldername, ntasks, 'greg.freeman@utexas.edu' );
events.updateTask= @updateTaskHtmlPost;

events.setup_command='cd  ~/experiment/natural_optim/;setup_natural_optim;';


testRunner(paramset,events);

% Experiment
function [outputImage,results]=getExperiment(image,settings) 


% precompute lasso image reconstruction
sz=size(image);

params=struct();
params.qmf =MakeONFilter('Symmlet',8);
params.reconstruct = 'lasso_tfocs';
params.save_wavelet_data=1;
params.save_coder_data=1;      
params.delta=settings.delta;
lambda1=0.001;
lambda2=0.1;
nIterationsOuter=15;

[intermediateImage,results]=CoarseFineExperiment(image,params);

W_op=results.W_op;
samples=results.samples;
vfine=results.vfine;
Phi=results.coder.Phi;
Npixels=numel(intermediateImage);
Ncoarse=Npixels-results.Nfine;
start_i=Ncoarse+1;
fine_samples = samples(start_i:end);

levels=5;
nBands=3*levels+1;
% startBand=8;
startBand=settings.startBand;

% load prior
prior_data=load('l1l2_pristine_stats');
if settings.cheat_r
    prior=prior_data.stats{settings.image};
else
    prior=prior_data.stats_wo_image{settings.image};   
end

eval= weighted_l1_l2_evaluator(Phi,fine_samples,vfine);

x_0=W_op*intermediateImage(:);
x=x_0;

options=struct();
options.f=eval.f;
options.df=eval.df;
options.x_0=x;
options.threshold=1e-2;
options.delta=settings.delta;
options.exact_linesearch=0;
options.max_iter=20;
options.threshold=1e-4;
options.precompute_image=false;
options.plot_linesearch=0;
options.plot_bt_linesearch=0;
options.save_gradient_components=0;
options.save_linesearch=0;
options.save_bt_linesearch=0;
options.bt_linesearch=1;
options.alpha = 0.3;

% save convergence
results.l1_converge=zeros(nBands,nIterationsOuter);
results.l2_converge=zeros(nBands,nIterationsOuter);

results.lambda1_band=zeros(nBands,nIterationsOuter);
results.lambda2_band=zeros(nBands,nIterationsOuter);
results.fx=zeros(1,nIterationsOuter);
results.f_components=zeros(3,nIterationsOuter);

% initialize lambdas
lambda1_band=ones(nBands,1)*lambda1;
lambda2_band=ones(nBands,1)*lambda2;
lambda1_cell=cell(nBands,1);
lambda2_cell=cell(nBands,1);

img=reshape(x,sz);
cellx=wavelet_band_tree(img,levels);

l1_goal=zeros(nBands,1);
l2_goal=zeros(nBands,1);
% estimate target l1 based on r
for iBand=startBand:nBands
    n=numel(cellx{iBand});
    l2_goal(iBand)=prior(iBand).l2;
    l1_goal(iBand)=prior(iBand).l2*sqrt(prior(iBand).r*n);
    %l1_goal(iBand)=prior(iBand).l1;
end

% initialize lambda vectors
for iBand=1:startBand-1
    lambda1_cell{iBand}=ones(size(cellx{iBand}))*lambda1;
    lambda2_cell{iBand}=zeros(size(cellx{iBand}));
end

% normalize with gsm model
if strcmp(settings.cheat_z,'yes')
    x_true=W_op*image(:);
    z_true=1./normalize_image_wavelet(x_true,levels,sz);
    z=z_true;
elseif strcmp(settings.cheat_z,'no')
    z=1./normalize_image_wavelet(x,levels,sz);    
end

for iIteration=1:nIterationsOuter  % outer problem       
    % target stats: weight l1 l2 relative to target r=f(alpha)
    if strcmp(settings.cheat_z,'em')
        % reestimate z every iteration
        z=1./normalize_image_wavelet(x,levels,sz);    
    end
    u2d=reshape(x.*z,sz);
    cellu=wavelet_band_tree(u2d,levels);
    l1_test=zeros(nBands,1);
    l2_test=zeros(nBands,1);
    for iBand=startBand:nBands
        l1_test(iBand)=norm(cellu{iBand}(:),1);
        l2_test(iBand)=norm(cellu{iBand}(:),2);
        
        lambda1_band(iBand) = lambda1_band(iBand)*(l1_test(iBand)/l1_goal(iBand))^2;
        lambda2_band(iBand) = lambda2_band(iBand)*(l2_test(iBand)/l2_goal(iBand))^2;
        
        lambda1_cell{iBand}=lambda1_band(iBand)*ones(size(cellx{iBand}));
        lambda2_cell{iBand}=lambda2_band(iBand)*ones(size(cellx{iBand}));
    end
    results.l1_converge(:,iIteration)=l1_test;
    results.l2_converge(:,iIteration)=l2_test;
    results.lambda1_band(:,iIteration)=lambda1_band;
    results.lambda2_band(:,iIteration)=lambda2_band;

    display(iIteration)
    lambda1=wavelet_band_image(lambda1_cell,levels);
    lambda2=wavelet_band_image(lambda2_cell,levels);
    lambda1=lambda1(:);
    lambda2=lambda2(:);

    % solve inner problem
    eval.setLambda(lambda1.*z,lambda2.*z);  % adjust weighting for evaluator
    [fs,ts,xs]=gradient_decent(options);
    x=xs(:,end);
    options.x_0=x; % set starting point for next iteration
    
    outputImage=reshape(W_op'*x,sz);
    
    figure(1),imagesc(outputImage)
    colormap gray
    figure(2),imagesc(intermediateImage)
    colormap gray
    
    [results.fx(iIteration),f_components]=eval.f(x);
    results.f_components(:,iIteration)=cell2mat(f_components)';

end
[results.fx_0,results.f_components_0]=eval.f(x_0);

results.intermediateImage=intermediateImage;
results.l2_goal=l2_goal;
results.l1_goal=l1_goal;

