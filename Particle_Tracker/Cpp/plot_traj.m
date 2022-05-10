clear;

% trajData = importdata('traj.txt','\t');
sliceData = importdata('slice_negtime.txt','\t');

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
minVal = -0.005e-2;
maxVal = 0.005e-2;
% maxVal = max(max(sliceData(sliceData < maxVal)));
% minVal = min(min(sliceData(sliceData > minVal)));
[minVal, maxVal] = bounds(sliceData(sliceData < maxVal & sliceData > minVal),'all');
sliceHeatMap = zeros(size(sliceData,1), numBins);
for i = 1:size(sliceData,1)
    h = histogram(sliceData(i,:), linspace(minVal, maxVal, numBins+1));
    sliceHeatMap(i,:) = h.Values;
end

figure(3)
% clim = [0 400];
imagesc(linspace(-10,100,size(sliceHeatMap,1)),...
    linspace(minVal, maxVal, h.NumBins),...
    sliceHeatMap')%, clim)
set(gca, 'Ydir', 'normal')
xlabel('t (ps)')
ylabel('Y (m)')
title('p_{y} = \pm 2 eV')
c = colorbar;
c.Label.String = 'Electron Count';
set(gca,'FontSize',20)

sliceSum = sum(sliceHeatMap, 2);
figure(4)
plot(linspace(-10,100,size(sliceHeatMap,1)), sliceSum);
xlim([-10, 100])
xlabel('t (ps)')
ylabel('Electron Count')
