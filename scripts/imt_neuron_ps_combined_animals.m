% imt_neuron_ps_combined_animals.m
function analysis_periods = imt_neuron_ps_combined_animals(structure,monkeys,task)
    analysis_periods = {};
    for mi=1:length(monkeys)
        analysis_periods{mi} = imt_neuron_ps(structure,monkeys(mi),task);
    end
end