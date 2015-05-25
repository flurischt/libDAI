figure;
hold on; grid on
x = logspace(-3,0);


bytesCycle = 4;
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

%loglog(opsCycle, intensity, 'ob');
%text(opsCycle,intensity,{' b003 baseline'},'VerticalAlignment','bottom');
loglog(intensity, opsCycle, 'ob');
text(intensity,opsCycle,{' b003 baseline'},'VerticalAlignment','bottom');

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

%loglog(opsCycle, intensity, 'ob');
%text(opsCycle,intensity,{' b003 baseline'},'VerticalAlignment','bottom');
loglog(intensity, opsCycle, 'ob');
text(intensity,opsCycle,{' b030 float'},'VerticalAlignment','bottom');
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

%loglog(opsCycle, intensity, 'ob');
%text(opsCycle,intensity,{' b003 baseline'},'VerticalAlignment','bottom');
loglog(intensity, opsCycle, 'ob');
text(intensity,opsCycle,{' b030 double'},'VerticalAlignment','bottom');
%% b032 single precision repeated
cycles = 54856082284;
% #total flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
% FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + FP_COMP_OPS_EXE.X87
flops = 1500002250 + 642000963 + 690001035;
runtime = 16.661; % seconds
bytes_sec = 1.071 * 1000000000;
bytes = runtime * bytes_sec;

bandwidthCycle = bytes/cycles;
opsCycle = flops/cycles;
intensity = opsCycle/bandwidthCycle;

%loglog(opsCycle, intensity, 'ob');
%text(opsCycle,intensity,{' b003 baseline'},'VerticalAlignment','bottom');
loglog(intensity, opsCycle, 'ob');
text(intensity,opsCycle,{' b032 float'},'VerticalAlignment','bottom');
%% b032 single precision repeated with DAI_RECOMMENDER_BOOST
cycles = 5796008694;
% #total flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
% FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + FP_COMP_OPS_EXE.X87
flops = 1494002241 + 570000855 + 702001053;
runtime = 16.234; % seconds
bytes_sec = 1.100 * 1000000000;
bytes = runtime * bytes_sec;

bandwidthCycle = bytes/cycles;
opsCycle = flops/cycles;
intensity = opsCycle/bandwidthCycle;

%loglog(opsCycle, intensity, 'ob');
%text(opsCycle,intensity,{' b003 baseline'},'VerticalAlignment','bottom');
loglog(intensity, opsCycle, 'ob');
text(intensity,opsCycle,{' b032 float BOOST'},'VerticalAlignment','bottom');
%% b032 double precision
cycles = 63334095001;
% #total flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
% FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + FP_COMP_OPS_EXE.X87
flops = 942001413 + 0 + 318000477;
runtime = 19.236; % seconds
bytes_sec = 5.328 * 1000000000;
bytes = runtime * bytes_sec;

bandwidthCycle = bytes/cycles;
opsCycle = flops/cycles;
intensity = opsCycle/bandwidthCycle;

%loglog(opsCycle, intensity, 'ob');
%text(opsCycle,intensity,{' b003 baseline'},'VerticalAlignment','bottom');
loglog(intensity, opsCycle, 'ob');
text(intensity,opsCycle,{' b032 double'},'VerticalAlignment','bottom');
%%
set(gca,'XScale','log')
set(gca,'YScale','log')
legend('Peak Performance', 'Maximal Memory Bandwidth','Location','northwest')
title('Roofline Plot')
xlabel('Operational Intensity [flops/byte]')
ylabel('Performance [flops/cycle]')
print(gcf, '-r150', 'roofline_single_vs_double.png', '-dpng');