%% All measurements are for B003 (Baseline)
data_points = [];

%% Performance Analysis 100k full
num_cycles = 377602017046; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 537756;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 14364021546 + 0 + 2664003996;
runtime = 114.685; % seconds
bandwidth_sec = 1.156; % GB/s

flops_cycle = num_flops/num_cycles;
data_points = [data_points; 1638, flops_cycle];

%% Performance Analysis u1
num_cycles = 351542527313; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 442075;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 13794020691 + 0 + 2364003546;
runtime = 106.770; % seconds
bandwidth_sec = 1.411; % GB/s

flops_cycle = num_flops/num_cycles;
% do not include in plot since it's basically the same as 100k
% data_points = [data_points; 1683, flops_cycle];

%% Performance Analysis [dataset user 405] 734 movies, user 405
num_cycles = 234140351210; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 330576;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 8928013392 + 0 + 1500002250;
runtime = 71.113; % seconds
bandwidth_sec = 1.153; % GB/s

flops_cycle = num_flops/num_cycles;
data_points = [data_points; 734, flops_cycle];

%% Performance Analysis [dataset 512] 272 movies, user 1-10
num_cycles = 1553744330613; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 2000409;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 66078099117 + 0 + 10242015363;
runtime = 471.622; % seconds
bandwidth_sec = 0.988; % GB/s

flops_cycle = num_flops/num_cycles;
data_points = [data_points; 272, flops_cycle];

%% Performance Analysis 32 movies, user 1-10
num_cycles = 72650108975; % CPU_CLK_UNHALTED.REF_TSC
num_messages = 163805;
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
num_flops = 3708005562 + 0 + 378000567;
runtime = 22.065; % seconds
bandwidth_sec = 0.567; % GB/s

flops_cycle = num_flops/num_cycles;
data_points = [data_points; 32, flops_cycle];
%% Plotting
plot(data_points(:,1), data_points(:,2), 'LineWidth',2);
hold on; grid on;

fontsize = 14;
%title('Baseline Performance for different Numbers of Movies', 'fontsize', fontsize);
xlabel('Number of Movies in Dataset', 'fontsize', fontsize)
y = ylabel('Performance [flops/cycle]', 'fontsize', fontsize,'rot', 0);
set(y, 'position', [500,0.0583,0]);
set(gca,'XTick',[32, 272, 734, 1638]);
set(gca,'fontsize', fontsize);
print(gcf, '-r150', 'performance_baseline.png', '-dpng');
saveas(gcf,'perf_baseline', 'pdf');
