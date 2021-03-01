kernels = {'original', 'x_EBC', 'x_RA', 'x_EBC_RA'};

load('../../nuPSOv16.mat', 'matrices', 'comm', 'info');
comm_real = comm; clear comm;
for k = 1:length(kernels)
    if exist(['../../results/nuPSOv16_' kernels{k} '_MCL_NMI.mat'],'file')
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
                            comm = mcl_for_graphs_v7(matrices{i1,i2,i3,i4,i5,i6}, C, kernels{k}, ['nPSO_MCL_' num2str(i6)], 'C:\cygwin64\bin');
                            
                            % NMI evaluation
                            if length(comm_real{i1,i2,i3,i4,i5,i6})/C < 100
                                NMI(i1,i2,i3,i4,i5,i6) = normalized_mutual_information(comm_real{i1,i2,i3,i4,i5,i6}, comm, 'adjusted');
                            else
                                NMI(i1,i2,i3,i4,i5,i6) = normalized_mutual_information(comm_real{i1,i2,i3,i4,i5,i6}, comm, 'unadjusted');
                            end
                        end
                        fprintf('%d %d %d %d %d [%.0fs]\n', i1, i2, i3, i4, i5, round(toc(time)));
                    end
                end
            end
        end
    end
    save(['../../results/nuPSOv16_' kernels{k} '_MCL_NMI.mat'], 'NMI')
end