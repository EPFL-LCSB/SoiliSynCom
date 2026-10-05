function updatedGEM = AddRavenGeneAnnotations(fastaFileName, GEM)

    % Read the fasta file and store annotations in a map
    geneAnnotations = containers.Map;
    fid = fopen(fastaFileName, 'r');
    tline = fgets(fid);
    while ischar(tline)
        if startsWith(tline, '>')
            geneId = extractBetween(tline, '>', ' ');
            geneName = extractAfter(tline, ' ');
            geneAnnotations(geneId{1}) = strtrim(geneName);
        end
        tline = fgets(fid);
    end
    fclose(fid);

    % Initialize a new field in GEM to store original gene loci
    GEM.geneLoci = GEM.genes;

    % Update the gene names in the GEM using the annotations map
    for i = 1:length(GEM.genes)
        originalGeneId = GEM.genes{i};
        if isKey(geneAnnotations, originalGeneId)
            GEM.genes{i} = geneAnnotations(originalGeneId);
        else
            warning('Annotation for gene %s not found', originalGeneId);
        end
    end
    
    % Assign the updated GEM to the output
    updatedGEM = GEM;
end
