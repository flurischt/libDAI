%% All measurements are for B003 (Baseline)
roofline_points = [];

%% Performance Analysis 100k full
num_cycles = 377602017046; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 537756;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 14364021546 + 0 + 2664003996;
runtime = 114.685; % seconds
bandwidth_sec = 1.156; % GB/s

bandwidth_cycle = (runtime*bandwidth_sec*1000000000)/num_cycles;
flops_cycle = num_flops/num_cycles;
intensity = flops_cycle/bandwidth_cycle;
roofline_points = [roofline_points; intensity, flops_cycle];

%% Performance Analysis [dataset user 405] 734 movies, user 405
num_cycles = 234140351210; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 330576;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 8928013392 + 0 + 1500002250;
runtime = 71.113; % seconds
bandwidth_sec = 1.153; % GB/s

bandwidth_cycle = (runtime*bandwidth_sec*1000000000)/num_cycles;
flops_cycle = num_flops/num_cycles;
intensity = flops_cycle/bandwidth_cycle;
roofline_points = [roofline_points; intensity, flops_cycle];

%% Performance Analysis [dataset 512] 272 movies, user 1-10
num_cycles = 1553744330613; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 2000409;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 66078099117 + 0 + 10242015363;
runtime = 471.622; % seconds
bandwidth_sec = 0.988; % GB/s

bandwidth_cycle = (runtime*bandwidth_sec*1000000000)/num_cycles;
flops_cycle = num_flops/num_cycles;
intensity = flops_cycle/bandwidth_cycle;
roofline_points = [roofline_points; intensity, flops_cycle];

%% Performance Analysis 32 movies, user 1-10
num_cycles = 72650108975; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 163805;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 3708005562 + 0 + 378000567;
runtime = 22.065; % seconds
bandwidth_sec = 0.567; % GB/s

bandwidth_cycle = (runtime*bandwidth_sec*1000000000)/num_cycles;
flops_cycle = num_flops/num_cycles;
intensity = flops_cycle/bandwidth_cycle;
roofline_points = [roofline_points; intensity, flops_cycle];

%% Roofline Plot Baseline
figure;
hold on; grid on
fontsize = 14;
x = logspace(-2,0);

maxBytesCycle = 4;
% flops / cycle
peakPerf = 2;
oneOpCycle = 1;

% Scalar peak performance
ys = ones(length(x),1).*peakPerf;
perfBound = loglog(x,ys,'-k','LineWidth',2);

% memory bandwidth
yb = x.*maxBytesCycle;
memBound = loglog(x,yb,'-k','LineWidth',2);
set(memBound(1),'color',[0.8 0.1 0.2]);
% measurements
p=loglog(roofline_points(:,1), roofline_points(:,2), '-o', 'LineWidth',2);
set(p(1),'color',[0.2 0.2 0.6])


legend('Peak Performance','Maximal Memory Bandwidth','Location','northwest');
xlabel('Operational Intensity [flops/byte]', 'fontsize', fontsize)
y=ylabel({'Performance [flops/cycle]'}, 'fontsize', fontsize,'rot', 0);
set(y, 'position', [0.0225,10.2,0]);
set(gca,'XScale','log')
set(gca,'YScale','log')
% t=title('Roofline Plot Baseline','fontsize', fontsize+1);
% titlepos = get(t,'position') + [0, 3, 0];
% set(t, 'position', titlepos);
set(gca,'fontsize', fontsize);
print(gcf, '-r150', 'roofline_baseline.png', '-dpng');

%% Measurements for B040
% compiled with gcc 4.5.2, -O3, cygwin64
roofline_points = [];
%% Performance Analysis 100k full, single precision
num_cycles = 3356005034; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 54912;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 72000108 + 12000018 + 18000027;
runtime = 1.019; % seconds
bandwidth_sec = 1.229; % GB/s

bandwidth_cycle = (runtime*bandwidth_sec*1000000000)/num_cycles;
flops_cycle = num_flops/num_cycles;
intensity = flops_cycle/bandwidth_cycle;
roofline_points = [roofline_points; intensity, flops_cycle];
%% Performance Analysis 100k full, double precision
num_cycles = 31454047181; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 448130;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 816001224 + 0 + 384000576;
runtime = 9.553; % seconds
bandwidth_sec = 2.176; % GB/s

bandwidth_cycle = (runtime*bandwidth_sec*1000000000)/num_cycles;
flops_cycle = num_flops/num_cycles;
intensity = flops_cycle/bandwidth_cycle;
roofline_points = [roofline_points; intensity, flops_cycle];
%% Performance Analysis [dataset user 405] 734 movies, user 405, single precision
num_cycles = 2538003807; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 54330;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 78000117 + 6000009 + 18000027;
runtime = 0.771; % seconds
bandwidth_sec = 1.568; % GB/s

bandwidth_cycle = (runtime*bandwidth_sec*1000000000)/num_cycles;
flops_cycle = num_flops/num_cycles;
intensity = flops_cycle/bandwidth_cycle;
roofline_points = [roofline_points; intensity, flops_cycle];
%% Performance Analysis [dataset user 405] 734 movies, user 405, double precision
num_cycles = 16308024462; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 275480;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 438000657 + 0 + 222000333;
runtime = 4.953; % seconds
bandwidth_sec = 2.238; % GB/s

bandwidth_cycle = (runtime*bandwidth_sec*1000000000)/num_cycles;
flops_cycle = num_flops/num_cycles;
intensity = flops_cycle/bandwidth_cycle;
roofline_points = [roofline_points; intensity, flops_cycle];
%% Performance Analysis [dataset 512] 272 movies, user 1-10, single precision
num_cycles = 9630014445; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 105670;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 204000306 + 30000045 + 102000153;
runtime = 2.925; % seconds
bandwidth_sec = 1.699; % GB/s

bandwidth_cycle = (runtime*bandwidth_sec*1000000000)/num_cycles;
flops_cycle = num_flops/num_cycles;
intensity = flops_cycle/bandwidth_cycle;
roofline_points = [roofline_points; intensity, flops_cycle];
%% Performance Analysis [dataset 512] 272 movies, user 1-10, double precision
num_cycles = 93388140082; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 2083777;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 2814004221 + 0 + 1356002034;
runtime = 30.253; % seconds
bandwidth_sec = 2.237; % GB/s

bandwidth_cycle = (runtime*bandwidth_sec*1000000000)/num_cycles;
flops_cycle = num_flops/num_cycles;
intensity = flops_cycle/bandwidth_cycle;
roofline_points = [roofline_points; intensity, flops_cycle];
%% Roofline Plot Baseline & B040
% Peak division perf
ys = ones(length(x),1).*(1/22);
divBound = loglog(x,ys,'--k','LineWidth',2);
% title('Roofline Plot Build 40','fontsize', fontsize);
legend([perfBound, divBound, memBound],{'Peak Performance', 'Peak Division Performance', 'Maximal Memory Bandwidth'},'Location','northwest');
% measurements
p=loglog(roofline_points(:,1), roofline_points(:,2), '-o', 'LineWidth',2);
set(p(1),'color',[0.1 0.5 0.1]);
text(0.2, 0.025,{'Baseline'},'VerticalAlignment','middle','HorizontalAlignment','center', 'fontsize', fontsize);
text(0.075, 0.025,{'Build 40'},'VerticalAlignment','middle','HorizontalAlignment','center', 'fontsize', fontsize);
print(gcf, '-r150', 'roofline_mixed.png', '-dpng');
saveas(gcf, 'roofline_mixed', 'pdf')