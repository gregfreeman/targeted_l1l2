function targeted_l1l2_one

param1.field='delta';
param1.values=num2cell(0.2);

param2.field='image';
param2.values=num2cell(1);

param3.field='cheat_r';
param3.values={false};

param4.field='cheat_z';
param4.values={'em'};


paramset=[param1,param2,param3,param4];

events=struct();
events.runExperiment=@targeted_l1l2;
events.loadInputData=@(settings) CSEvaluationImageGetter(settings);
events.evaluateMetrics= @extendedQualityMetrics;
events.storeOutputData= @pngStoreOutputData;
% events.startTasks = @( foldername, ntasks ) startTasksHtmlPost( foldername, ntasks, 'greg.freeman@utexas.edu' );
% events.updateTask= @updateTaskHtmlPost;

events.setup_command='cd  ~/experiment/natural_optim/;setup_natural_optim;';


testRunner(paramset,events);

