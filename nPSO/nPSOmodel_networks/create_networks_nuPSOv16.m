filename = 'nuPSOv16.mat';

N = [100,500,1000];
m = 2:2:16;
T = 0.1:0.2:0.7;
C = 3:3:12;
gamma = [2.01,2.5,3];
iters = 10;

matrices = cell(length(N),length(m),length(T),length(C),length(gamma),iters);
coords = cell(length(N),length(m),length(T),length(C),length(gamma),iters);
comm = cell(length(N),length(m),length(T),length(C),length(gamma),iters);

info = sprintf(['[D1] N = [100,500,1000]\n',...
    '[D2] m = 2:2:16\n',...
    '[D3] T = 0.1:0.2:0.7\n',...
    '[D4] C = 3:3:12\n',...
    '[D5] gamma = [2.01,2.5,3]\n',...
    'iters = 10']);

for i1 = 1:length(N)
    for i2 = 1:length(m)
        for i3 = 1:length(T)
            for i4 = 1:length(C)
                for i5 = 1:length(gamma)
                    display(num2str([i1 i2 i3 i4 i5]))
                    parfor i6 = 1:iters
                        [x, x_coords, x_comm] = nPSO_model(N(i1), m(i2), T(i3), gamma(i5), C(i4));
                        matrices{i1,i2,i3,i4,i5,i6} = x;
                        coords{i1,i2,i3,i4,i5,i6} = x_coords;
                        comm{i1,i2,i3,i4,i5,i6} = x_comm;
                    end
                end
            end
        end
    end
end

save(filename, 'matrices', 'coords', 'comm', 'info', '-v7.3')
