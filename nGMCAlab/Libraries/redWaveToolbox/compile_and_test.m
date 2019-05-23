mexingRoutine(['sources/' redWave])
mexingRoutine('redWaveThresholding')
mexingRoutine('redWaveInversion')

% add paths for tests
currentDir = pwd();
path([currentDir '/bin'], path)

% convenient function
al = @(x) x(:);

% redWave function must be used as:
% y = redWave(x, forward, options);
% with: - x: the data on which to apply the transfom (can be up to 4
%       dimensional data, this limitation can be modified in the c++ code)
%       - forward: computes the forward transform if >0, backaward
%       otherwise
%       - options: structure which can take the following values defined in
%       the next section. default value are provided with the following
%       error message:
fprintf(1, 'Test of error management (error only trigger a text because of\n an incompatibility of Matlab errors with win64):\n')
redWave()
%% wavelets options
clear options
%should not have any impact since all the fields are filled
options.verbose = 1; 
% number of scales
options.L = 3; 
% isometric mode for the wavelet (rice wavelet are not isometric, the
% inverse is not directly the transpose, but needs reweighting)
options.isometric = 0; 
% dimension on which to compute the wavelets
options.dimensions = 2; 
% wavelet filter. Filters need to have specific properties. Use
% MakeONFilter of the Wavelab toolbox for more choice.
options.qmf = [-.075765714789341 -.029635527645954 .497618667632458 .803738751805216 .297857795605542 -.099219543576935 -.012603967262261 .032223100604071];

%% example of use on 2 signals at once
S = zeros(2, 512);
S(1, 130) = 1; S(1, 170) = 0.4; S(1, 350) = 0.9;
S(2, 80) = 1; S(2, 200) = 0.4; S(2, 450) = 0.9;
kernel = exp( - abs(linspace( - 7, 7, 128))); 
S(1, :) = conv(S(1, :), kernel, 'same');
S(2, :) = conv(S(2, :), kernel, 'same');
plot(redWave(S, 1, options)', 'LineWidth', 2)
axis([0, 512 * 4, 0, 3]);
set(gca, 'Box', 'on', 'TickDir', 'out', 'XColor', 'black', 'YColor', 'black', 'YTick', [], 'XTick', [])

    
%% 1D comparison with wavelets from Rice University
fprintf('\n1D comparison with wavelets from Rice University\n');
if exist('mrdwt', 'file') == 3
    x = rand(16, 512);
    tic
    for k = 1 : 100
        y1 = redWave(x, 1, options);
    end
    t1 = toc();
    tic
    for k = 1 : 100
        y2 = riceWave(x, 1, options);
    end
    t2 = toc();
    
    fprintf(1, 'Forward acceleration (1D): x%1.2f.\n', t2 / t1);
    % must be of the order of the machine precision
    if norm(y1 - y2, Inf)<10^-9
        fprintf('1D forward transform tested and correct.\n');
    else
        error('Incorrect 1D forward transform.\n')
    end
    
    
    x1 = redWave(y1, -1, options);
    x2 = riceWave(y1, -1, options);
    
    if max(abs(al(x1 - x2))) < 10^-9
        fprintf('1D backward transform tested and correct.\n');
    else
        error('Incorrect 1D backward transform.\n')
    end
    
    
    
    %%
    
    %% 2D wavelets transform
    fprintf(1, '2D comparison with wavelets from Rice University.\n')
    
    clear options
    options.L = 3;
    options.isometric = 0;
    options.verbose = 1; %should not have any impact since all the fields are filled
    options.dimensions = [1 2];
    remaining_dim = 3;
    options.qmf = [-.075765714789341 -.029635527645954 .497618667632458 .803738751805216 .297857795605542 -.099219543576935 -.012603967262261 .032223100604071];
    x = rand(128, 256, 4);
    
    
    tic
    for k = 1 : 10
        y1 = riceWave(x, 1, options);
    end
    t2 = toc();
    
    
    tic
    for k = 1 : 10
        y2 = redWave(x, 1, options);
    end
    t1 = toc();
    y2rice = reshape2Dwavelets(y2, options);
    
    fprintf(1, 'Forward acceleration (2D): x%1.2f.\n', t2 / t1);
    
    % must be of the order of the machine precision
    if norm(al(y1 - y2rice), Inf)<10^ - 9
        fprintf('2D forward transform tested and correct.\n');
    else
        error('Incorrect 2D forward transform.\n')
    end
    
    
    x1 = riceWave(y1, -1, options);
    y1red = reshape2Dwavelets(y1, options, - 1);
    x2 = redWave(y1red, - 1, options);
    
    
    % must be of the order of the machine precision
    if max(abs(al(x1 - x2))) < 10^-9
        fprintf('2D backward transform tested and correct.\n');
    else
        error('Incorrect 2D backward transform.\n')
    end
else
    fprintf('Skipped checking against Rice wavelets (not found in the current path).\n')
end

%% test isometric mode
x = rand(128, 256, 1);
options.isometric = 1;
options.L = 3;
options.dimensions = [1 2];
y = redWave(x, 1, options);
yrice = reshape2Dwavelets(y, options); %careful: function not fully implemented
% in isometric mode, norme of y and x should be the same
% be careful however, for simplicity's sake, y has redundant
% coarse scales which are not used for the reconstruction
% and which must not be counted in the norm (hence the 
% reshaping).
if abs(norm(al(yrice)) - norm(al(x)))/ norm(al(x)) < 10^-9
    fprintf('Isometric mode tested and correct.\n');
else
    error('Isometric mode not working properly.\n')
end

%% example of 2D wavelets visualization
load example_lena
options.dimensions = [2, 1]; % easier to visualize in this order
imagesc(redWave(lena, 1, options));
colormap(linspace(0, 1, 256)' * [1, 1, 1]); set(gca,'DataAspectRatio',[1 1 1]);
set(gca, 'Box', 'on', 'TickDir', 'out', 'XColor', 'black', 'YColor', 'black', 'YTick', [], 'XTick', [])
% notice the 3 coarse versions of Lena
% only the top left one is used for the reconstruction


%%
%% Analysis proximal operator
fprintf(1, '\nAnalysis proximal operator\n');
fprintf(1, 'Pure C++ code Vs C++ / Matlab code.\n');
% computes the proximal operator of ||lambda .* (S W^T)||_1
% using the Forward-Backward algorithm
% details of parameters are provided in the matlab version:
% redWaveThresholding_m

clear options
options.verbose = 1; %should not have any impact since all the fields are filled
options.qmf = [-.075765714789341 -.029635527645954 .497618667632458 .803738751805216 .297857795605542 -.099219543576935 -.012603967262261 .032223100604071];
options.L = 1;
options.isometric = 1;
options.nonNegative = 1; %can handle both regular redundant wavelet thresholding and non-negative redundant wavelet thresholding
options.MaximumIteration = 20;
options.dimensions = [2, 3];
x = rand(4, 256, 256);
wsize = size(x);
wsize(options.dimensions) = wsize(options.dimensions) * (options.L + 1);
lambda = 5 * rand(wsize);


tic
for k = 1 : 4
    y1 = redWaveThresholding_m(x, lambda, options);
end
t1 = toc();

tic
for k = 1 : 4
    y2 = redWaveThresholding(x, lambda, options);
end
t2 = toc();


fprintf(1, 'C++ code acceleration: x%1.2f.\n', t1/t2);
% must be of the order of the machine precision
if max(abs(al(y1 - y2))) < 10^-9
    fprintf('Analysis proximal operator tested and correct.\n');
else
    error('Incorrect analysis proximal operator.\n')
end





%%
%% Inversion with non-negative and wavelet sparse priors
fprintf(1, '\nInversion with non-negative and wavelet sparse priors\n');
fprintf(1, 'Pure C++ code Vs C++ / Matlab code.\n');
% Solves argmin_{S>=0} ||Y-A * S||_2^2 + ||lambda .* (S W^T)||_1
% using Chambolle-Pock algorithm
% details of parameters are provided in the matlab version:
% redWaveNnInversion_m

% define the wavelets
clear options
options.verbose = 0; %should not have any impact since all the fields are filled
options.L = 3;
options.dimensions = 2;
options.MaximumIteration = 20;
options.isometric = 1;
options.nonNegative = 1; %can handle both regular redundant wavelet thresholding and non-negative redundant wavelet thresholding
options.qmf = [-.075765714789341 -.029635527645954 .497618667632458 .803738751805216 .297857795605542 -.099219543576935 -.012603967262261 .032223100604071];
W = @(x) redWave(x, 1, options);
Wt = @(x) redWave(x, -1, options);


% generate data
m = 32;
r = 8;
n = 512;
reference.A = rand(m, r);
reference.S = rand(r, n) .* (rand(r, n) > 0.99);
kernel = exp( - abs(linspace( - 7, 7, 128))); plot(kernel)
temp = reference.S';
for k = 1 : r;
    reference.S(k, :) = conv(reference.S(k, :), kernel, 'same');
end
plot(reference.S')
reference.Y = reference.A * reference.S;
reference.N = 0.01 * randn(m, n);

% parameters

Y = reference.A * (reference.S - 0.2) + reference.N;
S0 = 0 * reference.S;
A = reference.A;
lambdas = 0.2 * ones(size(W(S0)));
cost = @(x) 0.5 * norm(al(Y - A * x))^2 + sum(abs(al(lambdas .* (W(x)))));



tic
for k = 1 : 20
    Sc = redWaveInversion(S0, A' * A, A' * Y, lambdas, options);
end
t1 = toc();


tic
for k = 1 : 20
    Sm = redWaveNnInversion_m(S0, A' * A, A' * Y, lambdas, options);
end
t2 = toc();
plot(Sm')

fprintf(1, 'C++ code acceleration: x%1.2f.\n', t2 / t1);
% must be of the order of the machine precision
if max(abs(al(Sc - Sm))) < 10^-9
    fprintf('Analysis non-negative inversion tested and correct.\n');
else
    error('Incorrect analysis non - negative inversion.\n')
end


%%
%%
fprintf('\nThe toolbox is compiled and correct.\n')
fprintf('From now on, only add the bin folder in your path.\n')
fprintf('A list of wavelet filters which can be used for these wavelets.\n')
fprintf('can be found in\nthe MakeONFilter function of the WaveLab toolbox.\n')


