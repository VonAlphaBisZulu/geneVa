% Genome-scale
load iJO1366.mat model
cnap = CNAcobra2cna(model);
cnap = CNAsetGenericReactionData_with_array(cnap,'subSystems',model.subSystems);
cnap = CNAsetGenericReactionData_with_array(cnap,'geneProductAssociation', model.grRules);
% Remove this reaction to line up indices with exchange reaction list in iJO1366_extras
cnap = CNAdeleteReaction(cnap,find(ismember(cnap.reacID,{'Ec_biomass_iJO1366_WT_53p95M'})));

cnap.mue = find(ismember(cnap.reacID,{'Ec_biomass_iJO1366_core_53p95M'}));

grRules = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
load('iJO1366GeneNames.mat');
for i = ecoliGeneNames'
    grRules = strrep(grRules,i(1),i(2));
end

% A list of all products producible from glucose in iJO1366
load target_reactions
target_reactions = rxns_passed_preprocessing;

target = 'EX_etoh(e)'

del_exchanges = [];
if strmatch('iJO1366',model.description)
    load iJO1366_extras reac_off del_exchanges   
    % remove the product from deleted exchanges
    del_exchanges = setdiff(del_exchanges,find(ismember(cnap.reacID,{target})));
end

timeout = 600;	% s

solution_type = 3; % Search for feasible solutions. Optimal solution is returned if infinite timeout provided

% This defines the type of algo to run, can explicitely define a number of
% valves, or max number of valves (unweighted), or choose a number of
% valves based on a weight. Should probably use weighted first to see if
% any solution exists
objective_type = 'unweighted';  % 'weighted' or 'unweighted'

    % if weighted
    objective_weight = 1;

    % if unweighted
    number_of_valves = 3;
    valve_constraint = 'less_than_or_equal'; % 'equal' or 'less_than_or_equal'

growth_stage_min_mue = 0.5;	% percent of max growth rate

% use a low production_stage_min_yield when checking 
% whether coupling is possible at all
production_stage_min_yield = 0.5; % percent of max growth rate
production_stage_min_mue = 0.1; % biomass yield (gdw/mol)

% Number of threads to use, 0 for auto
num_threads = 0;

% Plot the production envelope of the result
% Note no visible envelope will be displayed if mu = 0
plot_envelope = true;

% Array of reaction indices which should not be used as valves. Empty string for none
no_valve_genes = {'s0001', 'spontanous'};

% notknockable genes
no_KO_genes = {'s0001', 'spontanous'};

% Alteratively a single_valve id can be chosen, only this valve will be allowed
single_valve = '';  % e.g. 'EX_o2(e)'

% reactions that can be knocked out
KO_reacs = {'EX_o2(e)'};
valve_reacs = {'EX_o2(e)'};

% Calculate and return solution (this is solved via JVM)
[kos, valves] = calc_geneVa(cnap, target, 'EX_glc(e)', growth_stage_min_mue, ...
    production_stage_min_yield, production_stage_min_mue, objective_type, ...
    objective_weight, number_of_valves, valve_constraint, solution_type, ...
    timeout, num_threads, single_valve, no_valve_genes, valve_reacs, plot_envelope, false,...
    reac_off, del_exchanges, no_KO_genes, KO_reacs, grRules);
