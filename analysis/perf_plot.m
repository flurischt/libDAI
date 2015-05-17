peak_perf = 2;

data_points = zeros(6,2);
data_points2 = zeros(6,2);

N = 32;
instr = 620944236;
cycles = 512333536;
flops = 24000060;
perf = flops/cycles;
ops_cycle = instr/cycles;
data_points2(1,:) = [N, ops_cycle];
data_points(1,:) = [N, perf];

N = 64;
instr = 1576639914;
cycles = 1354432651;
flops = 48000072;
perf = flops/cycles;
ops_cycle = instr/cycles;
data_points2(2,:) = [N, ops_cycle];
data_points(2,:) = [N, perf];

N = 128;
instr = 5800032921;
cycles = 6033692795;
flops = 240000648;
perf = flops/cycles;
ops_cycle = instr/cycles;
data_points2(3,:) = [N, ops_cycle];
data_points(3,:) = [N, perf];

N = 256;
instr = 10856396983;
cycles = 15026119602;
flops = 480001200;
perf = flops/cycles;
ops_cycle = instr/cycles;
data_points2(4,:) = [N, ops_cycle];
data_points(4,:) = [N, perf];

N = 512;
instr = 13481328085;
cycles = 18020959931;
flops = 600001404;
perf = flops/cycles;
ops_cycle = instr/cycles;
data_points2(5,:) = [N, ops_cycle];
data_points(5,:) = [N, perf];

N = 1024;
instr = 13477232228;
cycles = 17860077573;
flops = 600001212;
perf = flops/cycles;
ops_cycle = instr/cycles;
data_points2(6,:) = [N, ops_cycle];
data_points(6,:) = [N, perf];

plot(data_points(:,1), data_points(:,2));
hold on; grid on;

title('Performance Plot for different Numbers of Movies');
xlabel('Number of Movies in Dataset')
ylabel('Performance [flops/cycle]')

set(gca,'XTick',[32 64 128 256 512 1024]);
print(gcf, '-r150', 'performance_flops.png', '-dpng');


figure;
plot(data_points2(:,1), data_points2(:,2));
hold on; grid on;

title('Operations/Cycle Plot for different Number of Movies');
xlabel('Number of Movies in Dataset')
ylabel('[ops/cycle]')

set(gca,'XTick',[32 64 128 256 512 1024]);
print(gcf, '-r150', 'performance_ops.png', '-dpng');