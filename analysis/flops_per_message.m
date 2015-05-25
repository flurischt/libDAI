data = [];
data_points = 0;
builds = {};

%%
% number of cycles = CPU_CLK_UNHALTED.REF_TSC
% number of flops = FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE +
%                   FP_COMP_OPS_EXE.SSE_SCALAR_SINGLE + 
%                   FP_COMP_OPS_EXE.X87
%% B003 [Baseline] doubles, N=1
% -O3 -g
num_messages = 442075;
num_flops = 13974020961 + 0 + 2334003501;
num_cycles = 355138532707; % reported by vTune
num_cycles_print = 358095519066; % reported by the program itself
runtime_sec = 107.862;
bandwidth_GB_sec = 1.597;

data = [data;data_points, num_flops/num_messages];
data_points = data_points + 1;
builds = [builds, {'Baseline'}];

%% B030 floats, N=10
% -O3 -g
num_messages = 884068;
num_flops = 2052003078 + 576000864 + 816001224;
num_cycles = 59584089376; % reported by vTune
num_cycles_print = 5873471722; % reported by the program itself
runtime_sec = 18.097;
bandwidth_GB_sec = 0.931;

data = [data;data_points, num_flops/num_messages];
data_points = data_points + 1;
builds = [builds, {'B030'}];

%% B032 floats, N=10
% -O3 -g
num_messages = 884068;
num_flops = 1470002205 + 582000873 + 714001071;
num_cycles = 52960079440; % reported by vTune
num_cycles_print = 5584797774; % reported by the program itself
runtime_sec = 16.085;
bandwidth_GB_sec = 0.941;

data = [data;data_points, num_flops/num_messages];
data_points = data_points + 1;
builds = [builds, {'B032'}];

%% B036 floats, N=10
% -O3 -g
num_messages = 260215;
num_flops = 456000684 + 60000090 + 210000315;
num_cycles = 23580035370;
runtime_sec = 7.162;
bandwidth_GB_sec = 1.169;

data = [data;data_points, num_flops/num_messages];
data_points = data_points + 1;
builds = [builds, {'B036'}];

%% B042 floats, N=10
% -O3 -g
num_messages = 260215;
num_flops = 420000630 + 66000099 + 222000333;
num_cycles = 21540032310;
runtime_sec = 6.542;
bandwidth_GB_sec = 1.293;

data = [data;data_points, num_flops/num_messages];
data_points = data_points + 1;
builds = [builds, {'B040'}];

%% Plot styling
figure; hold on; grid on;
plot(data(:,1), data(:,2));
title('Flops per Message');
xlabel('Build')
ylabel('[flops/message]')
set(gca,'XTick',0:1:(data_points-1));
set(gca,'XTickLabel',builds);
print(gcf, '-r150', 'flops_per_message.png', '-dpng');