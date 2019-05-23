%% parameters
colors={'cyan', 'blue', 'red', 'green'};
algos={'oracle', 'nGMCA^S', 'nGMCA^H', 'nGMCA^{naive}'};
fontSize = 13;


%% load benchmark
clear options;
load bench_m_p9MC12_11Nov2013_1809_IEEE15dBa30;


%% draw options
options.algorithmNames = algos;
options.mode = 0;
options.colorList = colors;
options.fontSize = fontSize;
options.criteriaNames = {'SDR_S'};


%% draw
drawBenchmark(benchmark, options);
set(gca, ...
    'YTick'       , 0 : 4 : 24, ...
    'XTick'       , 0 : 40 : 200 );
axis([7, 210, 0, 22.5]);
