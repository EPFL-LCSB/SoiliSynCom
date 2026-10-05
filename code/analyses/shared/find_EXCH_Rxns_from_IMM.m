function EXCH_Rxns_alt = find_EXCH_Rxns_from_IMM(model_after_imm,DPs,num_min_nutrients) 
    
    ind_bfuse_M = getAllVar(model_after_imm,{'BFUSE'});
 
    %Extracting the BFUSE matrices of interest
    BFUSE_all = DPs(ind_bfuse_M,:);
    BFUSE_null = (BFUSE_all < 0.5);
    EXCH_Rxns_alt = cell(length(find(BFUSE_null(:,1))),length(BFUSE_null(1,:)));
 
    % Determining the names of the 31 open EXCHANGE rxns for each of the 1792 alternative
    % media
    for i=1:length(BFUSE_null(1,:)) 
        if length(find(BFUSE_null(:,i))) ~= num_min_nutrients
            error("Minimal medium is not composed of expected minimum number of nutrients, contrary to what was expected")
        else
            ind_rxns = ind_bfuse_M(BFUSE_null(:,i));
            EXCH_Rxns_alt(:,i) = model_after_imm.varNames(ind_rxns);
            for j=1:length(EXCH_Rxns_alt(:,i))
                EXCH_Rxns_alt{j,i} = strrep(EXCH_Rxns_alt{j,i},'BFUSE_','');
            end
        end
    end
    
end

