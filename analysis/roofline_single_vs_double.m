figure;
hold on; grid on
x = logspace(-2,0);


bytesCycle = 4.2;
% flops / cycle
peakPerf = 2;
oneOpCycle = 1;

% Scalar peak performance
ys = ones(length(x),1).*peakPerf;
loglog(x,ys,'-b');

% memory bandwidth
yb = x.*bytesCycle;
loglog(x,yb,'r');

%% Baseline
% b003
cycles = 339904509856;
% #total flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
% FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + FP_COMP_OPS_EXE.X87
flops = 14094021141 + 0 + 2304003456;
runtime = 103.235; % seconds
bytes_sec = 0.786 * 1000000000;
bytes = runtime * bytes_sec;


bandwidthCycle = bytes/cycles;
opsCycle = flops/cycles;
intensity = opsCycle/bandwidthCycle;

loglog(opsCycle, intensity, 'ob');
text(opsCycle,intensity,{' b003 baseline'},'VerticalAlignment','bottom');

%% b030 single precision
cycles = 6100009150;
% #total flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
% FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + FP_COMP_OPS_EXE.X87
flops = 186000279 + 48000072 + 60000090;
runtime = 1.853; % seconds
bytes_sec = 0.816 * 1000000000;
bytes = runtime * bytes_sec;

bandwidthCycle = bytes/cycles;
opsCycle = flops/cycles;
intensity = opsCycle/bandwidthCycle;

loglog(opsCycle, intensity, 'ob');
text(opsCycle,intensity,{' b030 single precision'},'VerticalAlignment','top');

%% b030 double precision
cycles = 62842094263;
% #total flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
% FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + FP_COMP_OPS_EXE.X87
flops = 1494002241 + 0 + 426000639;
runtime = 19.086; % seconds
bytes_sec = 2.082 * 1000000000;
bytes = runtime * bytes_sec;

bandwidthCycle = bytes/cycles;
opsCycle = flops/cycles;
intensity = opsCycle/bandwidthCycle;

loglog(opsCycle, intensity, 'ob');
text(opsCycle,intensity,{' b030 double precision'},'VerticalAlignment','top');

%%
set(gca,'XScale','log')
set(gca,'YScale','log')
legend('Peak Performance', 'Maximal Memory Bandwidth (fitted, 4.2 bytes/cycle)','Location','northwest')
title('Roofline Plot')
xlabel('Operational Intensity [flops/byte]')
ylabel('Performance [flops/cycle]')
print(gcf, '-r150', 'roofline_single_vs_double.png', '-dpng');