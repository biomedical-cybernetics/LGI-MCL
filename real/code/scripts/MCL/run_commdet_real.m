networks = {'karate', 'opsahl_8', 'opsahl_9', 'opsahl_10', 'opsahl_11', 'polbooks', 'football', 'polblogs'};

kernels = {'original','x_EBC_RA','x_RA','x_EBC_RA'};

for i = 1:length(networks)
    display(networks{i})
    load(['../../../networks/' networks{i} '.mat'], 'x', 'comm_real');
    C = numel(unique(comm_real)); % number of communities
    
    NMI_kernels = NaN(length(kernels),1);
    parfor k = 1:length(kernels)
        if exist(['../../results/' networks{i} '_' kernels{k} '_MCL_NMI.mat'],'file')
            continue
        end
        
        % community detection
        comm = mcl_for_graphs_v7(x, C, kernels{k}, ['real_MCL_' num2str(k)], 'C:\cygwin64\bin');
        
        % NMI evaluation
        if length(comm_real)/C < 100
            NMI_kernels(k) = normalized_mutual_information(comm_real, comm, 'adjusted');
        else
            NMI_kernels(k) = normalized_mutual_information(comm_real, comm, 'unadjusted');
        end
    end
    for k = 1:length(kernels)
        NMI = NMI_kernels(k);
        save(['../../results/' networks{i} '_' kernels{k} '_MCL_NMI.mat'], 'NMI')
    end
end