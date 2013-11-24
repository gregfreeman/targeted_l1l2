function targeted_l1l2_varybands_small

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
events.runExperiment=@targeted_l1l2;
events.loadInputData=@(settings) CSEvaluationImageGetter(settings);
events.evaluateMetrics= @extendedQualityMetrics;
events.storeOutputData= @pngStoreOutputData;
events.startTasks = @( foldername, ntasks ) startTasksHtmlPost( foldername, ntasks, 'greg.freeman@utexas.edu' );
events.updateTask= @updateTaskHtmlPost;

events.setup_command='cd  ~/experiment/natural_optim/;setup_natural_optim;';


testRunner(paramset,events);
