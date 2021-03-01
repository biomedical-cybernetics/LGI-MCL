kernels = {'original'};

load('../../nuPSOv16.mat', 'matrices', 'comm', 'info');
comm_real = comm; clear comm;
for k = 1:length(kernels)
    if exist(['../../results/nuPSOv16_' kernels{k} '_louvain_NMI.mat'],'file')
        continue
    end
    display(kernels{k})
    NMI = NaN(size(matrices));
    for i1 = 1:size(matrices,1)
        for i2 = 1:size(matrices,2)
            for i3 = 1:size(matrices,3)
                for i4 = 1:size(matrices,4)
                    for i5 = 1:size(matrices,5)
                        time = tic;
                        parfor i6 = 1:size(matrices,6)
                            C = length(unique(comm_real{i1,i2,i3,i4,i5,i6})); % number of communities
                            
                            % community detection
                            comm = louvain_R(matrices{i1,i2,i3,i4,i5,i6}, kernels{k}, ['nPSO_louvain_' num2str(i6)], 'C:\Users\mluser\Documents\R\R-3.4.3\bin\x64');
                            
                            % NMI evaluation
                            NMI_i6 = NaN(size(comm,2),1);
                            for l = 1:size(comm,2)
                                if length(comm_real{i1,i2,i3,i4,i5,i6})/C < 100
                                    NMI_i6(l) = normalized_mutual_information(comm_real{i1,i2,i3,i4,i5,i6}, comm(:,l), 'adjusted');
                                else
                                    NMI_i6(l) = normalized_mutual_information(comm_real{i1,i2,i3,i4,i5,i6}, comm(:,l), 'unadjusted');
                                end
                            end
                            NMI(i1,i2,i3,i4,i5,i6) = max(NMI_i6);
                        end
                        fprintf('%d %d %d %d %d [%.0fs]\n', i1, i2, i3, i4, i5, round(toc(time)));
                    end
                end
            end
        end
    end
    save(['../../results/nuPSOv16_' kernels{k} '_louvain_NMI.mat'], 'NMI')
end