%% All measurements are for B040
% compiled with gcc 4.5.2, -O3, cygwin64
data_points_single = [];
data_points_double = [];

%% Performance Analysis 100k full, single precision
num_cycles = 3356005034; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 54912;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 72000108 + 12000018 + 18000027;
runtime = 1.019; % seconds
bandwidth_sec = 1.229; % GB/s

flops_cycle = num_flops/num_cycles;
data_points_single = [data_points_single; 1638, flops_cycle];

%% Performance Analysis 100k full, double precision
num_cycles = 31454047181; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 448130;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 816001224 + 0 + 384000576;
runtime = 9.553; % seconds
bandwidth_sec = 2.176; % GB/s

flops_cycle = num_flops/num_cycles;
data_points_double = [data_points_double; 1638, flops_cycle];

%% Performance Analysis [dataset user 405] 734 movies, user 405, single precision
num_cycles = 2538003807; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 54330;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 78000117 + 6000009 + 18000027;
runtime = 0.771; % seconds
bandwidth_sec = 1.568; % GB/s

flops_cycle = num_flops/num_cycles;
data_points_single = [data_points_single; 734, flops_cycle];

%% Performance Analysis [dataset user 405] 734 movies, user 405, double precision
num_cycles = 16308024462; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 275480;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 438000657 + 0 + 222000333;
runtime = 4.953; % seconds
bandwidth_sec = 2.238; % GB/s

flops_cycle = num_flops/num_cycles;
data_points_double = [data_points_double; 734, flops_cycle];

%% Performance Analysis [dataset 512] 272 movies, user 1-10, single precision
num_cycles = 9630014445; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 105670;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 204000306 + 30000045 + 102000153;
runtime = 2.925; % seconds
bandwidth_sec = 1.699; % GB/s

flops_cycle = num_flops/num_cycles;
data_points_single = [data_points_single; 272, flops_cycle];

%% Performance Analysis [dataset 512] 272 movies, user 1-10, double precision
num_cycles = 93388140082; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 2083777;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 2814004221 + 0 + 1356002034;
runtime = 30.253; % seconds
bandwidth_sec = 2.237; % GB/s

flops_cycle = num_flops/num_cycles;
data_points_double = [data_points_double; 272, flops_cycle];

%% Performance Analysis 32 movies, user 1-100, single precision
num_cycles = 7810011715; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 46139;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 54000081 + 12000018 + 12000018;
runtime = 2.372; % seconds
bandwidth_sec = 0.935; % GB/s

flops_cycle = num_flops/num_cycles;
data_points_single = [data_points_single; 32, flops_cycle];

%% Performance Analysis 32 movies, user 1-100, double precision
num_cycles = 20396030594; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 872656;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 768001152 + 0 + 294000441;
runtime = 6.195; % seconds
bandwidth_sec = 0.04; % GB/s

flops_cycle = num_flops/num_cycles;
data_points_double = [data_points_double; 32, flops_cycle];
%% Plotting
plot(data_points_single(:,1), data_points_single(:,2),'-b', 'LineWidth',2);
hold on; grid on;
plot(data_points_double(:,1), data_points_double(:,2),'-r', 'LineWidth',2);

fontsize = 14;
title('B040 Performance for different Numbers of Movies', 'fontsize', fontsize);
legend('float','double', 'fontsize', fontsize);
xlabel('Number of Movies in Dataset', 'fontsize', fontsize)
ylabel('Performance [flops/cycle]', 'fontsize', fontsize)

set(gca,'XTick',[32, 272, 734, 1638]);
set(gca,'fontsize', fontsize);
print(gcf, '-r150', 'performance_b040.png', '-dpng');

