% imt_neuron_glm_combined_animals.m
function analysis_periodses = imt_neuron_glm_combined_animals(structure,monkeys,task,params)
    for mi=1:length(monkeys)
        analysis_periodses{mi} = imt_neuron_glm(structure,monkeys(mi),task,params);
    end
end