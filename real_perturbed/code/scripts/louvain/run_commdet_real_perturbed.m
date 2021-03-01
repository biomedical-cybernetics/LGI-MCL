networks = {'karate', 'opsahl_8', 'opsahl_9', 'opsahl_10', 'opsahl_11', 'polbooks', 'football', 'polblogs'};
perturb = {'linkrem10','linkadd10'};
kernels = {'original'};

for i = 1:length(networks)
    for p = 1:length(perturb)
        load(['../../../networks/' networks{i} '_' perturb{p} '.mat'], 'matrices', 'comm_real');
        C = numel(unique(comm_real)); % number of communities
        for k = 1:length(kernels)
            if exist(['../../results/' networks{i} '_' perturb{p} '_' kernels{k} '_louvain_NMI.mat'],'file')
                continue
            end
            fprintf(([networks{i} '_' perturb{p} ' - ' kernels{k}]));
            time = tic;
            NMI = NaN(length(matrices),1);
            parfor j = 1:length(matrices)

				% community detection
                comm = louvain_R(matrices{j}, kernels{k}, ['realpert_louvain_' num2str(j)], 'C:\Users\mluser\Documents\R\R-3.4.3\bin\x64');
                
                % NMI evaluation
                NMI_j = NaN(size(comm,2),1);
                for l = 1:size(comm,2)
                    if length(comm_real)/C < 100
                        NMI_j(l) = normalized_mutual_information(comm_real, comm(:,l), 'adjusted');
                    else
                        NMI_j(l) = normalized_mutual_information(comm_real, comm(:,l), 'unadjusted');
                    end
                end
                NMI(j) = max(NMI_j);
            end
            save(['../../results/' networks{i} '_' perturb{p} '_' kernels{k} '_louvain_NMI.mat'], 'NMI')
            time = round(toc(time));
            fprintf(' [%.0fs]\n', time)
        end
    end
end
