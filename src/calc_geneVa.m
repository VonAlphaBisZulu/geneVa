function [ kos, valves ] = calc_geneVa( cnap, target, glucose, growth_stage_min_mue, ...
    production_stage_min_yield, production_stage_min_mue, objective_type, ...
    objective_weight, number_of_valves, valve_constraint, solution_type, timeout,...
    num_threads, single_valve_reaction, no_valve_genes, valve_reacs, plot_envelope, distributed, ...
    reac_off, del_exchanges, no_KO_genes, KO_reacs, grRules)
    % MoVE Calculates multi-objective cutsets for a single product
    % Error codes
    % 0: Preparation error
    % 1: FVA error
    % 2: LP error
    % 3: Parameter value error
    % 4: No solution found from MILP
    % 9: Other error

    %% This is the configuration information for cplex's distributed search.
    %% These variables must be set for a distributed search
    if distributed
        node_1 = getenv('NODE_1')
        node_2 = getenv('NODE_2')
        node_3 = getenv('NODE_3')
        node_4 = getenv('NODE_4')

        vmconfig = sprintf(['<?xml version="1.0" encoding="US-ASCII"?>'...
                    '<vmc>'...
                    '  <machine name="machine1">'...
                    '    <transport type="TCP/IP">'...
                    '      <address host="%s" port="8050" />'...
                    '    </transport>'...
                    '  </machine>'...
                    '  <machine name="machine2">'...
                    '    <transport type="TCP/IP">'...
                    '      <address host="%s" port="8050" />'...
                    '    </transport>'...
                    '  </machine>'...
                    '  <machine name="machine3">'...
                    '    <transport type="TCP/IP">'...
                    '      <address host="%s" port="8050" />'...
                    '    </transport>'...
                    '  </machine>'...
                    '  <machine name="machine4">'...
                    '    <transport type="TCP/IP">'...
                    '      <address host="%s" port="8050" />'...
                    '    </transport>'...
                    '  </machine>'...
                    '</vmc>'],node_1,node_2,node_3,node_4);
    end

    total_solution = [];
    total_solution.solution = 0;
    total_solution.name = target;

    %% Set cplex parameters
    cplex_inner= setup_cplex_inner_class_access();
    emphasizeAccuracy= cplex_inner.ParameterSet.constructor.newInstance([]);
    emphasizeAccuracy.setParam(cplex_inner.BooleanParam.NumericalEmphasis, true);
    % emphasizeAccuracy.setParam(cplex_inner.BooleanParam.PreInd, false);
    emphasizeAccuracy.setParam(cplex_inner.DoubleParam.EpOpt, 1e-9);
    emphasizeAccuracy.setParam(cplex_inner.DoubleParam.EpRHS, 1e-9);
    % emphasizeAccuracy.setParam(cplex_inner.IntParam.ScaInd, 1);
    emphasizeAccuracy.setParam(cplex_inner.IntParam.AdvInd, 0);

    % explicit parameter set for solving the MILP
    MILPparameters= cplex_inner.ParameterSet.constructor.newInstance([]);
    MILPparameters.setParam(cplex_inner.IntParam.MIPEmphasis, cplex_inner.MIPEmphasis.Balanced);
    MILPparameters.setParam(cplex_inner.IntParam.MIPDisplay, 2);
    MILPparameters.setParam(cplex_inner.DoubleParam.EpInt, 0);
    MILPparameters.setParam(cplex_inner.DoubleParam.EpRHS, 1e-9);
    MILPparameters.setParam(cplex_inner.BooleanParam.NumericalEmphasis, true);
    MILPparameters.setParam(cplex_inner.DoubleParam.WorkMem, 1024*10);
    MILPparameters.setParam(cplex_inner.IntParam.RandomSeed, 2);
    MILPparameters.setParam(cplex_inner.IntParam.Threads,num_threads); % limit threads for testing
    MILPparameters.setParam(cplex_inner.IntParam.NodeFileInd, 3); % Saves nodefile to disk (2 for uncompressed, 3 for compressed)
    if distributed
        MILPparameters.setParam(cplex_inner.DistMIP.Rampup.TimeLimit, timeout/4);
    end
    num_evs_lim = 5; %max number of solutions

    %% Prepare the model
    mue = cnap.mue;

    single_valve_id = find(ismember(cnap.reacID,{single_valve_reaction}));
    disp(single_valve_id);
    glucose = find(ismember(cnap.reacID,{glucose}));
    atpm = find(ismember(cnap.reacID,{'ATPM'}));
    
    KO_reacs = ismember(cnap.reacID,KO_reacs);
    valve_reacs = ismember(cnap.reacID,valve_reacs);
     
    exchange_rxns = ~cellfun(@isempty,(regexp(cellstr(cnap.reacID),'EX_.*')));
    DM_rxns = ~cellfun(@isempty,(regexp(cellstr(cnap.reacID),'DM_.*')));
    transport_rxns = ~cellfun(@isempty,regexp(CNAgetGenericReactionData_as_array(cnap,'subSystems'),'Transport'));
    exchange_rxns(single_valve_id) = false;
    
    % remove exchange reactions
    cnap.reacMin([reac_off(:); del_exchanges(:)]) = 0;
    cnap.reacMax([reac_off(:); del_exchanges(:)]) = 0;

    product_index = find(ismember(cnap.reacID,{target}));
    product_spec_index = find(cnap.stoichMat(:,product_index));
    if length(product_spec_index) > 1
        total_solution.error = 0;
        total_solution.error_message = 'More than one product index found for the target exchange';
        disp('More than one product index found for the target exchange');
        return;
    end

    glucose_uptake_limit = -cnap.reacMin(glucose);
    min_atpm= cnap.reacMin(atpm);
    num_seeds= 1;

    %% Begin solution
    
    protected_reacs = DM_rxns | exchange_rxns | transport_rxns;
    protected_reacs(KO_reacs) = false;
    
    cnap.objFunc = zeros(cnap.numr, 1);
    cnap.objFunc(product_index) = -1;
    fv = CNAoptimizeFlux(cnap,[],[],2);
    max_prod= fv(product_index);
    gluc_up= -fv(glucose);
    
    cnap.objFunc = zeros(cnap.numr, 1);
    cnap.objFunc(mue) = -1;
    fv = CNAoptimizeFlux(cnap,[],[],2);
    max_mue= -fv(mue);    

    prod_min_yield= max_prod/gluc_up*production_stage_min_yield;
    min_mue_yield_prod_stage= -max_mue/gluc_up*production_stage_min_mue;
    
    %% Set D1, D2, T
    T= zeros(0, cnap.numr);
    T(1, glucose)= -1;
    T(2, [product_index, glucose])= [1, prod_min_yield]; % +prod_min_yield because glucose influx has a negative sign!
    T(3, atpm)= -1;
    t= [glucose_uptake_limit; 0; -min_atpm];
    % production stage
    D_pr_stage= zeros(0, cnap.numr);
    D_pr_stage(1, [mue, glucose])= [-1, -min_mue_yield_prod_stage]; % growth yield
    D_pr_stage(2, glucose)= -1;
    D_pr_stage(3, [product_index, glucose])= [-1, -prod_min_yield];
    D_pr_stage(4, atpm)= -1;
    d_pr_stage= [0; glucose_uptake_limit; 0; -min_atpm]; % growth yield
    % growth stage
    D_gr_stage= zeros(0, cnap.numr);
    D_gr_stage(1, mue)= -1;
    D_gr_stage(2, glucose)= -1;
    D_gr_stage(3, atpm)= -1;
    d_gr_stage= [max_mue*growth_stage_min_mue; glucose_uptake_limit; -min_atpm]; %valve
    
    if 1
    %% Gene Extension
        % get GR-rules in right format and remove notknockables
        [~, ~, genes, gpr_rules] = CNAgenerateGPRrules( cnap, grRules, 0);
        targetable_genes = ones(length(genes),1);
        targetable_genes(ismember(genes,no_valve_genes)) = 0; 
        targetable_genes(ismember(genes,no_KO_genes)) = 0;
        displ('FVA to determine blocked reactions.',1);
        [lb_fva , ub_fva ] = CNAfluxVariability(cnap,[],[],2);
        reac_off = find(lb_fva == 0 & ub_fva == 0);
        displ('FVA to determine essential reactions at production stage. Used for reducing gene protein reaction (GPR) rules.',1);
        [lb_Dp, ub_Dp] = CNAfluxVariability(cnap,[],[],2,1:cnap.numr,D_pr_stage,d_pr_stage);
        essential = sign((abs(lb_Dp)>cnap.epsilon).*lb_Dp) .* ... % when lower bound sign equals
                    sign((abs(ub_Dp)>cnap.epsilon).*ub_Dp) == 1  ;
        essential(protected_reacs) = true; % exclude also transporters etc.
        displ('Reducing set of GPR rules.',1);
        % remove rules for reactions that are off and continue with the reduced set
        gpr_rules = gpr_rules(~ismember([gpr_rules.reaction],reac_off));
        rule_removelist = zeros(1,length(gpr_rules));
        [reac_abund,reac] = hist([gpr_rules(:).reaction],unique([gpr_rules(:).reaction]));
        cgenes = {gpr_rules(:).genes}; % cell array of genes
        remaining_genes = intersect(cell2mat(cgenes),1:length(genes)); % all genes that weren't only catalyzing blocked reactions
        for i = remaining_genes
            % if a gene only catalyzes essential reactions -> set gene notknockable
            rules_i = cellfun(@(x) ismember(i,x),cgenes); % bool vector: Rules that contain the gene
            if all(essential([gpr_rules(rules_i).reaction]))
                targetable_genes(i) = 0;
            end
            % if the gene is essential to one essential reaction -> set gene notknockable
            for j = unique([gpr_rules(rules_i).reaction])
                                 % if all reaction-gene-rules for a particular reaction contain the gene
                if essential(j) && sum([gpr_rules(rules_i).reaction] == j) == reac_abund(j == reac)
                    targetable_genes(i) = 0;
                end
            end
            % if gene is notknockable remove it from all rules
            if targetable_genes(i) == 0 && any(rules_i)
                for j = find(rules_i)
                    gpr_rules(j).strGene = gpr_rules(j).strGene(gpr_rules(j).genes ~= i);
                    gpr_rules(j).genes = setdiff(gpr_rules(j).genes,i);
                end
            end
        end
        % if all genes of one rule are notknockable, delete all rules for the same reaction
        % because the reaction can never be knocked out.
        % Here, through preprocessing, some rules don't contain any genes anymore (genes = []) 
        % and are thus also notknockable
        for i = 1:length(gpr_rules)
            if all(targetable_genes(gpr_rules(i).genes) == 0) 
                rule_removelist([gpr_rules(:).reaction] == gpr_rules(i).reaction) = 1;
            end
        end
        gpr_rules = gpr_rules(~rule_removelist);
        % Remove non-minimal GPR rules
        gpr_rules = remove_nonminimal(gpr_rules);

        no_KO_genes = ismember(genes,no_KO_genes);
        no_valve_genes = ismember(genes,no_valve_genes);
        % genes that don't occur in rules anymore are notknockable
        targetable_genes(setdiff(1:length(genes),[gpr_rules(:).genes])) = 0;
        targetable_genes(targetable_genes==0) = nan;
        displ('Integrating GPR rules',1);
        [ cnap, rmap, gmap, ~, ~,~,~, rType ] = CNAintegrateGPRrules( cnap, gpr_rules, [],[],[],[],targetable_genes);
        D_pr_stage = D_pr_stage*rmap;
        D_gr_stage = D_gr_stage*rmap;
        T  =  T*rmap;

        no_KO_genes     = any(gmap(no_KO_genes,:),1);
        no_valve_genes  = any(gmap(no_valve_genes,:),1);
        KO_reacs        = any(rmap(KO_reacs,:),1);
        valve_reacs     = any(rmap(valve_reacs,:),1);
        reac_off        = any(rmap(reac_off,:),1);

        valvable = (rType == 'g');
        valvable(no_valve_genes) = 0;
        valvable(valve_reacs) = 1;
        knockable = (rType == 'g');
        knockable(no_KO_genes) = 0;
        knockable(KO_reacs) = 1;

        displ('Compressing network',1);
        non_compress_reacs = unique(find(any([D_pr_stage ; D_gr_stage ; T],1)))';
        [~,~,cmp_mapReac,~,cmp_cnap]=CNAcompressMFNetwork(cnap,non_compress_reacs,[],1,0,1,reac_off,0);
        cmp_D_gr_stage = D_gr_stage*cmp_mapReac;
        cmp_D_pr_stage = D_pr_stage*cmp_mapReac;
        cmp_T          = T *cmp_mapReac;
        cmp_knockable  = logical(abs(cmp_mapReac)'*knockable');
        cmp_valvable   = logical(abs(cmp_mapReac)'*valvable');

        %% bounds and knockables
        displ('FVA to bound target region',1);
        [ lb_fva , ub_fva ] = CNAfluxVariability(cmp_cnap,[],[],2);
        lb_relevant = abs(lb_fva-cmp_cnap.reacMin) <= cmp_cnap.epsilon;
        ub_relevant = abs(ub_fva-cmp_cnap.reacMax) <= cmp_cnap.epsilon;
        if any(lb_fva < 0 & lb_fva > -cmp_cnap.epsilon)
            warning('Target region: Some reactions might be irreversible but not treated as such. This has probably numerical reasons.');
        end
        T_entries = [-lb_relevant , ub_relevant]';
        T_bounds  = [-cmp_cnap.reacMin' ; cmp_cnap.reacMax']; % Using FVA bounds is identical: [-reacMin_FVA'; reacMax_FVA']; . Consider when numerical issues occur.
        [~,T_reac] = find(T_entries);
        cmp_T(size(cmp_T,1)+(1:length(T_reac)),:) = sparse(1:length(T_reac),T_reac,T_entries(T_entries~=0),length(T_reac),cmp_cnap.numr); % append target matrix
        cmp_t = t;
        cmp_t(    size(t,1) +   (1:length(T_reac)),:) = T_bounds(T_entries~=0); % append target vector

        displ('FVA to bound growth stage desired region',1);
        [lb_Dg, ub_Dg] = CNAfluxVariability(cmp_cnap,[],[],2,1:cmp_cnap.numr,cmp_D_gr_stage,d_gr_stage); 
        displ('FVA to bound production stage desired region',1);
        [lb_Dp, ub_Dp] = CNAfluxVariability(cmp_cnap,[],[],2,1:cmp_cnap.numr,cmp_D_pr_stage,d_pr_stage);

    else
        %% Reaction-based
        displ('FVA to find inactive reactions',1);
        [lb_fva, ub_fva] = CNAfluxVariability(cnap,[],[],2);
        non_compress_reacs = unique(find(any([D_pr_stage ; D_gr_stage ; T],1)))';
        rType = repmat('r',1,cnap.numr);
        rmap = eye(cnap.numr);
        gmap = zeros(cnap.numr,0);
        reac_off = find(lb_fva == 0 & ub_fva == 0);
        [~,~,cmp_mapReac,~,cmp_cnap]=CNAcompressMFNetwork(cnap,non_compress_reacs,[],1,0,1,reac_off,0);
        cmp_D_pr_stage = D_pr_stage*cmp_mapReac;
        cmp_D_gr_stage = D_gr_stage*cmp_mapReac;
        cmp_T  = T *cmp_mapReac;
        protected_reacs = ~logical(abs(cmp_mapReac)'*~protected_reacs);
        displ('FVA to bound growth stage desired region',1);
        [lb_Dg, ub_Dg] = CNAfluxVariability(cmp_cnap,[],[],2,1:cmp_cnap.numr,cmp_D_gr_stage,d_gr_stage); % relevant for static KOs
        displ('FVA to bound production stage desired region',1);
        [lb_Dp, ub_Dp] = CNAfluxVariability(cmp_cnap,[],[],2,1:cmp_cnap.numr,cmp_D_pr_stage,d_pr_stage);
        essential_prod_stage = (sign((abs(lb_Dp)>cnap.epsilon).*lb_Dp) .* ... % when lower bound sign equals
                                sign((abs(ub_Dp)>cnap.epsilon).*ub_Dp) == 1); % upper bound sign
        % define knockables and valves
        valvable = ones(cmp_cnap.numr,1);
        valvable(essential_prod_stage) = 0;
        valvable(protected_reacs) = 0; % exclude also transporters etc.

        knockable = ones(cmp_cnap.numr,1);
        knockable(protected_reacs) = 0; % exclude also transporters etc.
        knockable(~valvable) = 0;
    end
    r_irrev_pos = cmp_cnap.reacMin >= 0;
    
    % Check if knockables is a subset of valves
    if ~all(knockable<=valvable)
         error('Knockable candidates must be a subset of valve candidates');
    end
    
    %% Continue with MoVE
    
    s= 1;
    while s <= num_seeds %&& (~is_cs || ~is_des) && ~infeasible
        %     obj= ConstrainedMinimalCutSetsEnumerator(rd, irrev_rd, [], inh, ub, cuts, flux_lb{product_index}, flux_ub{product_index}, des, db);
        %     obj.cpx.setParam(cplex_inner.DoubleParam.EpGap, 0.98);
        obj= Orthogonal_cMCSEnumerator(cmp_cnap.stoichMat, r_irrev_pos', [], cmp_T, cmp_t, cmp_knockable,...
            lb_Dp, ub_Dp, cmp_D_pr_stage, d_pr_stage, lb_Dg, ub_Dg, cmp_D_gr_stage, d_gr_stage, find(~cmp_valvable));

        if distributed
            % obj.cpx.readCopyVMConfig(vmconfig); % set virtual machine configuration (or read from a file with obj.cpx.readVMConfig(vmcfile);)
            obj.cpx.copyVMConfig(vmconfig); % set virtual machine configuration (or read from a file with obj.cpx.readVMConfig(vmcfile);)
        end
        
        % Set-up the problem specifics
        if strcmp(objective_type,'unweighted')
            %     This is the case 1 where we want only one valve
            if strcmp(valve_constraint,'equal')
                % Equality constraint
                obj.cpx.addEq(obj.diff_expr, number_of_valves); %implements Sum(yold)-sum(ynew)=1
            elseif strcmp(valve_constraint,'less_than_or_equal')
                %Inequality constraint
                obj.cpx.addLe(obj.diff_expr, number_of_valves); %Implementing const. Sum(yold)-sum(ynew)<=1
            else
                total_solution.error = 3;
                total_solution.error_message = 'Unknown valve_constraint';
                error('Unknown valve_constraint');
            end

            % unweighted objective
            obj= set_objective_function(obj); %without weight
        elseif strcmp(objective_type,'weighted')
            %This is the case 2 where we can have more valves by adjusting weight
            obj= set_objective_function(obj,objective_weight);
        else
            total_solution.error = 3;
            total_solution.error_message = 'Unrecognized objective type';
            error('Unrecognized objective type');
        end

        % calculate cutsets
        % obj.cpx.exportModel('xyz.sav')
        [obj, mcs, status, ct]= shortest_minimal_cut_sets(obj, cmp_cnap.numr, num_evs_lim, timeout, solution_type, 0, MILPparameters);
        total_solution.status = status;
        if isempty(mcs)
            disp('no solution found');
            kos = [];
            valves = [];
            return
        end
        % %Following        
        proper_mcs= make_minimal_cs(cmp_cnap.stoichMat, cmp_cnap.reacMin >= 0, cmp_T, cmp_t, mcs, cmp_cnap.reacMin, cmp_cnap.reacMax);
        s= s + 1;
    end


    pause(1) % Ensure that this output comes after the Java stdout

    %% check if cut set is minimal (fulfilled when sum(mcs) == sum(proper_mcs))
    a = sum(mcs);
    b = sum(proper_mcs);
    if a > b
        disp('MCS not equal to proper MCS')
    end

    try
        fv2on= round(obj.cpx.getValues(obj.fv2_on));
    catch exception
        total_solution.error = 4;
        total_solution.error_message = 'Failed to calculate solution';
        disp('Failed to calculate solution')
        disp(target)
        disp(exception)
        disp(target);
        disp(exception.message);
        total_solution.name = target;
        total_solution.solution = 0;
        return
    end
    disp('Compressed Solution');
    disp('Valve: 1, KO:0');
    disp(cellstr(char(strcat(num2str(fv2on(mcs)), ':', cmp_cnap.reacID(mcs,:)))));
    % Print the generic solution (reduced model)
    kos = mcs & ~fv2on;
    valves = mcs & fv2on;
    disp('Verifying Valves and KOs');
    % checking growth stage
    if ~mcs_feas(cmp_cnap,kos,cmp_D_gr_stage,d_gr_stage)
        disp('Desired region for growth stage not feasible');
    end
    % checking production stage
    if mcs_feas(cmp_cnap,mcs,cmp_T,cmp_t) && ~mcs_feas(cmp_cnap,mcs,cmp_D_pr_stage,d_pr_stage)
        disp('Desired region for production stage not feasible or target region feasible');
    end
    
    % Expanding solution
    mcs_exp= expand_mcs(mcs, cmp_mapReac');
    valves_in_sol = cmp_mapReac*valves & valvable';
    valid_mcs = true(size(mcs_exp,2),1);
    disp(['The found (compressed) solution could be expanded to ' num2str(size(mcs_exp,2)) ' strain designs. Verifying strain designs.']);
    for i = 1:size(mcs_exp,2)
    % filter expanded MCS for knockables and valvables and number of valves
        if ~all(knockable(mcs_exp(:,i)) | valvable(mcs_exp(:,i)))
            valid_mcs(i) = 0;
            continue;
        end
        if sum(valves_in_sol(mcs_exp(:,i))) > number_of_valves
            disp('solution has more valves than allowed');
        end
    % filter expanded MCS for feasibility of growth and production stage
        kos = mcs_exp(:,i) & ~valves_in_sol;
        if ~mcs_feas(cnap,kos,D_gr_stage,d_gr_stage)
            valid_mcs(i) = 0;
            continue;
        end
    % checking properties in production stage
        if mcs_feas(cnap,mcs_exp(:,i),T,t) && ~mcs_feas(cnap,mcs_exp(:,i),D_pr_stage,d_pr_stage)
            valid_mcs(i) = 0;
            continue;
        end
    end
    mcs_exp = mcs_exp(:,find(valid_mcs,1));
    disp([num2str(size(mcs_exp,2)) ' valid strain designs.']);
    disp('Returning the knockouts and valves for the first solution:');
    
    %% First solution is prepared to be returned
    kos = mcs_exp & ~valves_in_sol;
    valves = mcs_exp & valves_in_sol;
    
    gene_kos = genes(logical((rType == 'g' & kos') * gmap'));
    if any(rType == 'r' & kos')
        reac_kos = cellstr(cnap.reacID(rType == 'r' & kos',:));
    else
        reac_kos = {};
    end
    gene_valves = genes(logical((rType == 'g' & valves') * gmap'));
    if any(rType == 'r' & valves')
        reac_valves = cellstr(cnap.reacID(rType == 'r' & valves',:));
    else
        reac_valves = {};
    end

    % Print the result
    disp('Knockouts');
    disp('  genes:');
    disp(gene_kos);
    disp('  reactions:');
    disp(reac_kos);
    disp('')
    disp('Valves');
    disp('  genes:');
    disp(gene_valves);
    disp('  reactions:');
    disp(reac_valves);
    
    if plot_envelope
        % Plot the three envelopes: wt, growth, production
        f = figure;
        hold on;
        
        mue = find(any(rmap(mue,:),1));
        glucose = find(any(rmap(glucose,:),1));

        % The substrate uptake must be fixed to ensure correct envelopes
        % Alternatively, we can calculate the yield explicitely
        cnap.reacMax(glucose) = cnap.reacMin(glucose);
        productionEnvelope(cnap,kos,'r',find(ismember(cnap.reacID,{target})),mue);
        productionEnvelope(cnap,kos | valves,'blue',find(ismember(cnap.reacID,{target})),mue);
        productionEnvelope(cnap,[],'k',find(ismember(cnap.reacID,{target})),mue);
        legend('growth stage','production stage','wild-type');
        title(strrep(target,'EX_',''));
    end
    pause(1);
end

function [xValues, targetValues, lineHandle] = productionEnvelope(cnap, deletions, lineColor, yreac, xreac)

cnap.reacMin(deletions) = 0;
cnap.reacMax(deletions) = 0;

if ~isfield(cnap,'macroDefault')
    cnap.macroDefault = 0;
end

% Run FBA to get upper bound for biomass
[solMin, solMax] = CNAfluxVariability(cnap,[],[],2,xreac);

% Create biomass range vector
xValues = linspace(solMin,solMax,10);

% Max/min for target production
fixed_fluxes = nan(cnap.numr,1);

for i = 1:length(xValues)
    fixed_fluxes(xreac) = xValues(i);
    [yLB(i), yUB(i)] = CNAfluxVariability(cnap,fixed_fluxes,[],2,yreac);
end

% Plot results
lineHandle=plot([xValues fliplr(xValues)],[yUB fliplr(yLB)],lineColor,'LineWidth',2);
axis tight;
%ylabel([strrep(targetRxn,'_','-') ' (mmol/gDW h)']);
%xlabel('Growth rate (1/h)');

xValues = xValues';
targetValues = [yLB' yUB'];
end

%% remove non-minimal rules
% A or (A and B) = A
function gpr_rules = remove_nonminimal(gpr_rules)
    maxreac = max([gpr_rules(:).reaction]);
    reac_rule = sparse([gpr_rules(:).reaction],1:length(gpr_rules),1,maxreac,length(gpr_rules));
    [~,~,rule_reac_abund] = unique(reac_rule','rows','stable'); 
    [rule_abund,isorule_map] = hist(rule_reac_abund,unique(rule_reac_abund));
    isorule_grule_ismin = false(1,length(gpr_rules));
    % 1. for non-isoenzymes, the only existing gene rule is the shortest
    isorule_grule_ismin( ismember(rule_reac_abund', isorule_map(rule_abund == 1)') ) = 1;
    % 2. for isoenzymes, gene rules are compared to eliminate redundancies
    for i = isorule_map(rule_abund > 1)'
        genes_i = {gpr_rules(rule_reac_abund == i).genes};
        isoenz = find(rule_reac_abund == i);
        isoenz_gene_rule = sparse(  cell2mat(genes_i),...
                                    repelem(1:length(isoenz), cellfun(@length,genes_i)),...
                                    1);
        % remove redundancies
        for j = 1:length(isoenz) % mark and eliminate redundant and non-minimal gene rules:
            if sum(all(repmat(isoenz_gene_rule(:,j),1,size(isoenz_gene_rule,2)) >= isoenz_gene_rule))==1  % A or (A and B) = A
                isorule_grule_ismin(isoenz(j)) = 1; % rule is minimal e.g. A
            else
                isorule_grule_ismin(isoenz(j)) = 0; % rule is not minimal e.g. (A and B)
            end
        end
    end
    gpr_rules = gpr_rules(isorule_grule_ismin); % remove redundant enzymes -> this also reduces set of genes
end

function feas = mcs_feas(cnap,mcs,V,v)
    cnap.reacMin(mcs) = 0;
    cnap.reacMax(mcs) = 0;
    feas = testRegionFeas(cnap,V,v,2);
end
