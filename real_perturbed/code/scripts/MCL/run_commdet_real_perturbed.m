networks = {'karate', 'opsahl_8', 'opsahl_9', 'opsahl_10', 'opsahl_11', 'polbooks', 'football', 'polblogs'};
perturb = {'linkrem10','linkadd10'};
kernels = {'original', 'x_EBC', 'x_RA', 'x_EBC_RA'};

for i = 1:length(networks)
    for p = 1:length(perturb)
        load(['../../../networks/' networks{i} '_' perturb{p} '.mat'], 'matrices', 'comm_real');
        C = numel(unique(comm_real)); % number of communities
        for k = 1:length(kernels)
            if exist(['../../results/' networks{i} '_' perturb{p} '_' kernels{k} '_MCL_NMI.mat'],'file')
                continue
            end
            fprintf(([networks{i} '_' perturb{p} ' - ' kernels{k}]));
            time = tic;
            NMI = NaN(length(matrices),1);
            parfor j = 1:length(matrices)

				% community detection
                comm = mcl_for_graphs_v7(matrices{j}, C, kernels{k}, ['realpert_MCL_' num2str(j)], 'C:\cygwin64\bin');
                
                % NMI evaluation
                if length(comm_real)/C < 100
                    NMI(j) = normalized_mutual_information(comm_real, comm, 'adjusted');
                else
                    NMI(j) = normalized_mutual_information(comm_real, comm, 'unadjusted');
                end
            end
            save(['../../results/' networks{i} '_' perturb{p} '_' kernels{k} '_MCL_NMI.mat'], 'NMI')
            time = round(toc(time));
            fprintf(' [%.0fs]\n', time)
        end
    end
end