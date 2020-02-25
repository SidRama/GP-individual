function pathway_rand_overlap(selected_genes, bf_score, save_file)

    %   Code to perform the permutation analysis and compute the null
    %   distribution of the adjusted geometric mean metric.
    %   Input:
    %       selected_genes - .mat file containing the lists of selected
    %                        gene expressions of all the case-control pairs 
    %                        in separate columns. The genes are sorted in 
    %                        the descending order of mean expression across
    %                        all case-control pairs.
    %       bf_score - chosen BF-score.
    %       file_save - resutls file name (ending with .mat)
    %   Output:
    %       .mat file containing the null distribution, i.e. the adjusted
    %       geometric mean computed for each of the 100000 label 
    %       permutations.

    children_mean_sorted = load(selected_genes);
    
    gm_stat = [];
    wt_stat = [];
    ck = 1;
    for set_id = 1:100000
        
        L = 19310;
        label_shuffle = randperm(L, L);
        labels = children_mean_sorted(:,1);
        labels = labels(label_shuffle, 1);
        children_mean_sorted(:,1) = labels;
        child_sel = cell(size(children_mean_sorted,2), 1)

        % obtain selected genes for each case-control pair
        for idx = 1:size(children_mean_sorted,2)
            child_sel(idx) = children_mean_sorted(find(cell2mat(children_mean_sorted(:, idx)) > bf_score), 1);
        end
        % finding the overlap per pathway
        pathway_count = 1;
        output = cell(16808, 10);
        for file_count = 1:17
            file_name = './data/msigdb_partition';
            file_name = strcat(file_name, num2str(file_count));
            file_name = strcat(file_name, '.mat');
            temp = load(file_name, 't');
            msigdb_partition =  temp.t;
            clear temp;
            for pathway_index = 1:size(msigdb_partition, 1)
                pathway = {};
                count = 0;
                 for idx = 3:size(msigdb_partition, 2)
                    if(strcmp(msigdb_partition{pathway_index, idx}, ''))
                        break;
                    else
                        count = count + 1;
                        pathway{count} = msigdb_partition{pathway_index, idx};
                    end      
                 end
                pathway = pathway';       
                % pathway comparison
                output{pathway_count, 1} = msigdb_partition{pathway_index,1};
                
                child_in = cell(size(children_mean_sorted,2), 1)
                for idx = 1:size(children_mean_sorted, 2)
                    child_in(idx) = intersect(child_sel(idx), pathway);
                    output{pathway_count, idx+1} = num2cell(length(child_in(idx)));
                end
                        
                pathway_count = pathway_count + 1;
            end
        end
 
        % get scores for null distribution
        % set alpha
        alpha = 0.5;
        
        % divide by number of enriched genes
        f_child = zeros(pathway_count-1, size(children_mean_sorted,2));
        for idx = 1:size(children_mean_sorted,2);
            f_child(:, idx) = (cell2mat(output(:, idx + 1)) + alpha)./length(child_sel(idx));
        
        mean_var = mean(var(f_child, 0, 2));
        
        for pathway_count = 1:size(output, 1)
            gm_stat(ck, pathway_count) = geomean(f_child(pathway_count, :));
            
        end
        ck = ck + 1;

    end

    save(save_file, 'gm_stat')
end
        
        
