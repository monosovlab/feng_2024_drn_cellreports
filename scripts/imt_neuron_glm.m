% imt_neuron_glm_offer.m
function analysis_periods = imt_neuron_glm(structure,monkey,task,params)
    savefilename = ['data//glm_offer_analysis_periods_' char(task) '_' char(structure) '_' char(monkey) '_params.offer.anawind' num2str(params.offer.anawind) '.mat'];
    load(savefilename);
end