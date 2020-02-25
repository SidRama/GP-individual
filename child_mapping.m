function child_mapping(probeset_ratio_file, annotation_file, probeset_names, save_file, genelistsort)
    
    %   Code to perform mapping from probe-sets to genes.
    %   Input:
    %      probeset_ratio_file - .mat file containing BF-scores computed from 
    %                            compute_ratios function.
    %      annotation_file - .mat file containing many-to-many mapping of
    %                        genes to probesets.
    %      probeset_names - .mat file containing probeset names.
    %      save_file - results file name (ending with .mat).
    %      genelistsort - .mat file containing a list of all genes in the 
    %                     pathway database in descending order of occurence
    %                     in the database
    %  Output:
    %      .mat file with N x 3 cell array. N is the number of genes. The
    %      first column contains the gene name, the second column contains
    %      the probe-set name and the third columns contains the corresponding
    %      BF-score. 

    genelistsort = load(genelistsort);
    annotation = load(annotation_file);
    ratios = load(probeset_ratio_file);
    ratios = ratios_noise_7;
    probe_name = load(probeset_names);
    count = 1;
    selected_probesets = {};
    for i = 1:size(ratios, 1)
        if (ratios{i} > 3)
            selected_probesets{count, 1} = probe_name{i, 1};
            selected_probesets{count, 2} = ratios{i};
            count = count + 1;
        end
    end

    % selected probesets obtained

    % map to gene labels
    gene_mapped = {};
    for cur_probe = 1:size(selected_probesets,1)
        flag = 0;
        for i = 1:size(annotation, 1)
            for j = 2:size(annotation, 2)
               if  (strcmp(annotation{i, j}, ''))
                   break;
               end
               if (strcmp(annotation{i, j},  selected_probesets{cur_probe, 1}))
                   gene_mapped{cur_probe, 1} = annotation{i, 1};
                   gene_mapped{cur_probe, 2} = selected_probesets{cur_probe, 1};
                   gene_mapped{cur_probe, 3} = selected_probesets{cur_probe, 2};
                   flag = 1;
                   break;
               end
            end
            if (flag == 1)
                break;
            end
        end
        if (mod(cur_probe, 100) == 0)
         display(strcat(num2str(cur_probe), ' probes mapped'));
        end
    end
    gene_map_processed = {};
    unique_gene_labels = unique(gene_mapped(:, 1));
    for idx = 1:size(unique_gene_labels, 1)
        max_ratio = -1;
        max_idx = -1;
        for gene_num = 1:size(gene_mapped,1)
            if(strcmp(gene_mapped{gene_num, 1}, unique_gene_labels{idx, 1}))
                if(gene_mapped{gene_num, 3} > max_ratio)
                    max_ratio = gene_mapped{gene_num, 3};
                    max_idx = gene_num;
                end
            end
        end
        gene_map_processed{idx, 1} = gene_mapped{max_idx, 1};
        gene_map_processed{idx, 2} = gene_mapped{max_idx, 2};
        gene_map_processed{idx, 3} = gene_mapped{max_idx, 3}; 
    end

    % comparing with pathway list
    gene_list_chosen = {};
    count = 0;
    not_found = 0;
    notfound_set = {};
    % Handle probe-set that maps to multiple gene-names
    for idx = 1:size(gene_map_processed,1)
           split = strsplit(gene_map_processed{idx, 1}, ' /// ');
            max_val = -1;
            max_idx = -1;
           if (size(split, 2) > 1)

                for i = 1:size(split, 2)
                    for g = 1:size(genelistsort, 1)
                        if(strcmp(split(1, i), genelistsort{g, 2}))
                            if( genelistsort{g, 1} > max_val)
                                max_val = genelistsort{g,1};
                                max_idx = g;
                            end
                            break;
                        end
                    end
                end
            if (max_val == -1)
                not_found = not_found + 1;
                notfound_set{not_found} = split;
            else
                count = count + 1;
                gene_list_chosen{count, 1} = genelistsort{max_idx, 2}; 
                gene_list_chosen{count, 2} = gene_map_processed{idx, 2};
                gene_list_chosen{count, 3} = gene_map_processed{idx, 3};
            end
           else

               for g = 1:size(genelistsort, 1)
                        if(strcmp(split(1, 1), genelistsort{g, 2}))
                            if( genelistsort{g, 1} > max_val)
                                max_val = genelistsort{g,1};
                                max_idx = g;
                            end
                            break;
                        end
               end

                if (max_val == -1)
                    not_found = not_found + 1;
                    notfound_set{not_found} = split;
                else
                    count = count + 1;
                    gene_list_chosen{count, 1} = genelistsort{max_idx, 2}; 
                    gene_list_chosen{count, 2} = gene_map_processed{idx, 2};
                    gene_list_chosen{count, 3} = gene_map_processed{idx, 3};
                end

           end
           if (mod(idx, 100) == 0)
               display(strcat(num2str(idx), ' processed'));
           end
    end

    % check for duplicates and retain only highest ratio
    unique_gene_labels = unique(gene_list_chosen(:, 1));
    child_list_chosen = {};
    for idx = 1:size(unique_gene_labels, 1)
        max_ratio = -1;
        max_idx = -1;
        for gene_num = 1:size(gene_list_chosen,1)
            if(strcmp(gene_list_chosen{gene_num, 1}, unique_gene_labels{idx, 1}))
                if(gene_list_chosen{gene_num, 3} > max_ratio)
                    max_ratio = gene_list_chosen{gene_num, 3};
                    max_idx = gene_num;
                end
            end
        end
        child_list_chosen{idx, 1} = gene_list_chosen{max_idx, 1};
        child_list_chosen{idx, 2} = gene_list_chosen{max_idx, 2};
        child_list_chosen{idx, 3} = gene_list_chosen{max_idx, 3}; 
    end

    save(save_file, 'child_list_chosen');

 end
