% imt neuron ps
function analysis_periods = imt_neuron_ps(structure,monkey,task)
  savefilename = ['data//ps_analysis_periods_'  char(task) '_' char(structure) '_' char(monkey) '.mat'];
  load(savefilename);
end