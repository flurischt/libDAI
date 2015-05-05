figure;
hold on; grid on
x = logspace(-1,1);


cyclesSec = 3300000000;
bytesSec = 7000000000;
%bytesCycle = bytesSec/cyclesSec;
bytesCycle = 4;
% ops / cycle
peakPerf = 4;
oneOpCycle = 1;

% Scalar peak performance
ys = ones(length(x),1).*peakPerf;
loglog(x,ys,'-b');

% memory bandwidth
yb = x.*bytesCycle;
loglog(x,yb,'r');

% fitted memory bandwidth
yb = x.*2.15;
loglog(x,yb,'-.r');

%% Baseline
% b003
cycles = (352756529134 + 352866529299)/2;
ops = (398286597429 + 398078597117)/2;
bytes = (174800000000 + 174500000000) / 2;

bandwidthCycle = bytes/cycles;
opsCycle = ops/cycles;
intensity = opsCycle/bandwidthCycle;

loglog(opsCycle, intensity, 'xk');
text(opsCycle,intensity,{' Baseline (b003)'});

%% Message Product Optimization
% b006
cycles = 241028361542;
ops = 186468279702;
bytes = 111500000000;

bandwidthCycle = bytes/cycles;
opsCycle = ops/cycles;
intensity = opsCycle/bandwidthCycle;

loglog(opsCycle, intensity, 'xk');
text(opsCycle,intensity,{' b006'});

%% b012
cycles = 227828341742;
ops = 161504242256;
bytes = 184900000000;

bandwidthCycle = bytes/cycles;
opsCycle = ops/cycles;
intensity = opsCycle/bandwidthCycle;

loglog(opsCycle, intensity, 'xk');
text(opsCycle,intensity,{' b012'});

%% b020
cycles = 82162123243;
ops = 28860043290;
bytes = 89480000000;

bandwidthCycle = bytes/cycles;
opsCycle = ops/cycles;
intensity = opsCycle/bandwidthCycle;

loglog(opsCycle, intensity, 'xk');
text(opsCycle,intensity,{' b020'},'VerticalAlignment','top');

%% b021
cycles = 80530120795;
ops = 28816043224;
bytes = 85960000000;

bandwidthCycle = bytes/cycles;
opsCycle = ops/cycles;
intensity = opsCycle/bandwidthCycle;

loglog(opsCycle, intensity, 'xk');
text(opsCycle,intensity,{' b021'}, 'VerticalAlignment','bottom');

%%
set(gca,'XScale','log')
set(gca,'YScale','log')
legend('Peak Performance', 'Maximal Memory Bandwidth', 'Fitted Memory Bandwidth (~2.1 bytes /cycle)','Location','northwest')
title('Roofline Plot')
xlabel('Operational Intensity [ops/byte]')
ylabel('Performance [ops/cycle]')
print(gcf, '-r150', 'roofline.png', '-dpng');