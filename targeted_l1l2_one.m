function targeted_l1l2_one

param1.field='delta';
param1.values=num2cell(0.2);

param2.field='image';
param2.values=num2cell(1);

paramset=[param1,param2];

events=struct();
events.runExperiment=@run_targeted_l1l2;
events.loadInputData=@(settings) CSEvaluationImageGetter(settings);
events.evaluateMetrics= @extendedQualityMetrics;
events.storeOutputData= @pngStoreOutputData;
% events.startTasks = @( foldername, ntasks ) startTasksHtmlPost( foldername, ntasks, 'greg.freeman@utexas.edu' );
% events.updateTask= @updateTaskHtmlPost;
events.setup_command='cd  ~/experiment/targeted_l1l2/;setup_targeted_l1l2;';

% load prior data
prior_data=load('l1l2_pristine_stats');

testRunner(paramset,events);

function [outputImage,results] = run_targeted_l1l2(image,settings) 
settings.prior_data=prior_data.stats_wo_image{settings.image};   
[outputImage,results] = targeted_l1l2(image,settings) ;
end

end
