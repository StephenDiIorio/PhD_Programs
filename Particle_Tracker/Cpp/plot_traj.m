clear;

% trajData = importdata('traj.txt','\t');
sliceData = importdata('slice_laminar.txt','\t');

n0 = 0.17863390738e26 * 1e-6;
n0_coeff = 531409.3265537234;
n0_const = n0_coeff / (100 * sqrt(n0));
% trajData = trajData * n0_const;
sliceData = sliceData * n0_const;

% figure(1);
% hold on;
% for i = 2 : 3 : size(trajData,2)
%     plot3(trajData(:,i),...
%         trajData(:,i+1),...
%         trajData(:,i+2));
% end
% hold off;
% xlabel('X (m)')
% ylabel('Y (m)')
% zlabel('Z (m)')
% set(gca,'FontSize',20)
%
% figure(2);
% final_plane = trajData(end, 3:3:size(trajData,2));
% histogram(final_plane, 100)
% xlabel('Y (m)')
% ylabel('Number of Hits')
% set(gca,'FontSize',20)

numBins = 100;
sliceHeatMap = zeros(size(sliceData,1), numBins);
for i = 1:size(sliceData,1)
    h = histogram(sliceData(i,:), numBins);
    sliceHeatMap(i,:) = h.Values;
end
[minVal, maxVal] = bounds(sliceData,'all');

figure(3)
clim = [0 200];
imagesc(linspace(0,100,size(sliceHeatMap,1)),...
    linspace(minVal, maxVal, h.NumBins),...
    sliceHeatMap', clim)
set(gca, 'Ydir', 'normal')
xlabel('t (ps)')
ylabel('Y (m)')
title('p_{y} = \pm 0 eV (Laminar Beam)')
c = colorbar;
c.Label.String = 'Electron Count';
set(gca,'FontSize',20)
