peak_perf = 2;

data_points = zeros(6,2);

N = 32;
instr = 620944236;
flops = 24000060;
perf = flops/instr;
data_points(1,:) = [N, perf];

N = 64;
instr = 1576639914;
flops = 48000072;
perf = flops/instr;
data_points(2,:) = [N, perf];

N = 128;
instr = 5800032921;
flops = 240000648;
perf = flops/instr;
data_points(3,:) = [N, perf];

N = 256;
instr = 10856396983;
flops = 480001200;
perf = flops/instr;
data_points(4,:) = [N, perf];

N = 512;
instr = 13481328085;
flops = 600001404;
perf = flops/instr;
data_points(5,:) = [N, perf];

N = 512;
instr = 13477232228;
flops = 600001212;
perf = flops/instr;
data_points(6,:) = [N, perf];

plot(data_points(:,1), data_points(:,2));
hold on;

title('Performance Plot for different Number of Movies');
xlabel('Number of Movies in Dataset')
ylabel('Performance [flops/cycle]')

set(gca,'XTick',[32 64 128 256 512 1024]);