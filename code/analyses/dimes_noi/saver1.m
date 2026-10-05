function r1=saver1(saving_directory_new,filename,model_p,DPsimm)
 % save workspace
            modeldescription=model_p.id;
            save(strcat(saving_directory_new,filename,'.mat'),'model_p','DPsimm');
            % save results
            filename = strcat(modeldescription,'_PhenoMappingExchanges');
            r1 = load(strcat(saving_directory_new,filename,'.mat'));
end