networks = {'karate', 'opsahl_8', 'opsahl_9', 'opsahl_10', 'opsahl_11', 'polbooks', 'football', 'polblogs'};
kernels = {'original'};

for i = 1:length(networks)
    display(networks{i})
    load(['../../../networks/' networks{i} '.mat'], 'x', 'comm_real');
    C = numel(unique(comm_real)); % number of communities
    
    NMI_kernels = NaN(length(kernels),1);
    for k = 1:length(kernels)
        if exist(['../../results/' networks{i} '_' kernels{k} '_infomap_NMI.mat'],'file')
            continue
        end
        
        % community detection
        comm = infomap(x, kernels{k}, ['real_infomap_' num2str(k)], 'C:\cygwin64\bin');
        
        % NMI evaluation
        NMI = NaN(size(comm,2),1);
        for l = 1:size(comm,2)
            if length(comm_real)/C < 100
                NMI(l) = normalized_mutual_information(comm_real, comm(:,l), 'adjusted');
            else
                NMI(l) = normalized_mutual_information(comm_real, comm(:,l), 'unadjusted');
            end
        end
        NMI_kernels(k) = max(NMI);
    end
    for k = 1:length(kernels)
        NMI = NMI_kernels(k);
        save(['../../results/' networks{i} '_' kernels{k} '_infomap_NMI.mat'], 'NMI')
    end
end