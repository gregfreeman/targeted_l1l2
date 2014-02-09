function targeted_l1l2_small

param1.field='delta';
% param1.values=num2cell(0.1:0.1:0.4);
param1.values=num2cell(0.05:0.05:0.45);

param2.field='image';
% param2.values=num2cell(3);
param2.values=num2cell(1:17);

param3.field='cheat_r';
% param3.values={false,true};
param3.values={false};

param4.field='cheat_z';
% param5.valutes={'no','yes','em'};
param4.values={'em'};

param5.field='image_options';
param5.values={{'cropsize',[128 128]}};

paramset=[param1,param2,param3,param4,param5];

events=struct();
events.runExperiment=@targeted_l1l2;
events.loadInputData=@(settings) CSEvaluationImageGetter(settings);
events.evaluateMetrics= @extendedQualityMetrics;
events.storeOutputData= @pngStoreOutputData;
events.startTasks = @( foldername, ntasks ) startTasksHtmlPost( foldername, ntasks, 'greg.freeman@utexas.edu' );
events.updateTask= @updateTaskHtmlPost;

events.setup_command='cd  ~/experiment/targeted_l1l2/;setup_targeted_l1l2;';


testRunner(paramset,events);

