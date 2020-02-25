function pathway_overlap(children_genes, bf_score, file_save)

    %   Code to compute the adjusted gemoetric mean score for each pathway
    %   using the true pathway overlaps.
    %   Input:
    %       children_genes - .mat file containing the list of selected
    %                                genes for each case-control pair.
    %       bf_score - chosen BF-score.
    %       file_save - results file name (ending with .mat) 
    %   Output:
    %       .mat file with adjusted geometric mean scores (enrichment scores)
    %       for each pathway.


    ck=1;
    pathway_count = 1;
    children_gene_selected = load(children_genes)
    output = cell(16808, size(children_gene_selected, 2));
    child_sel = cell(size(children_gene_selected,2), 1)

    % obtain selected genes for each case-control pair
    for idx = 1:size(children_gene_selected,2)
        child_sel(idx) = children_gene_selectedn_sorted(find(cell2mat(children_gene_selected(:, idx)) > bf_score), 1);
    end
    for file_count = 1:17
        fprintf('%d in progress\n',file_count);
        file_name = './data/msigdb_partition';
        file_name = strcat(file_name, num2str(file_count));
        file_name = strcat(file_name, '.mat');
        temp = load(file_name, 't');
        msigdb_partition =  temp.t;
        clear temp;
        for pathway_index = 1:size(msigdb_partition, 1)
            pathway = {};
            count = 0; % number of genes in pathway
             for idx = 3:size(msigdb_partition, 2)
                if(strcmp(msigdb_partition{pathway_index, idx}, ''))
                    break;
                else
                    count = count + 1;
                    pathway{count} = msigdb_partition{pathway_index, idx};
                end      
             end
            pathway = pathway';       
            % pathway
            output{pathway_count, 1} = msigdb_partition{pathway_index,1};
            child_in = cell(size(children_gene_selected,2), 1)
            for idx = 1:size(children_gene_selected, 2)
                child_in(idx) = intersect(child_sel(idx), pathway);
                output{pathway_count, idx+1} = num2cell(length(child_in(idx)));
            end
            pathway_count = pathway_count + 1;
        end
    end

    % get scores
    % set alpha
    alpha = 0.5;

    % divide by number of enriched genes
    f_child = zeros(pathway_count-1, size(children_gene_selected,2));
    for idx = 1:size(children_gene_selected,2);
        f_child(:, idx) = (cell2mat(output(:, idx + 1)) + alpha)./length(child_sel(idx));
    
    mean_var = mean(var(f_child, 0, 2));

    for pathway_count = 1:size(output, 1)
        gm_stat(ck, pathway_count) = geomean(f_child(pathway_count, :));
            
    end

    save(save_file, 'gm_stat')
