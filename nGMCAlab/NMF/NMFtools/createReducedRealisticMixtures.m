% createReducedRealisticMixtures.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jeremy Rapin (jeremy.rapin@cea.fr)
% Created on 20/08/13, last modified on 16/7/2014
% 
% This software is governed by the CeCILL  license under French law and
% abiding by the rules of distribution of free software.  You can  use, 
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info". 
% 
% As a counterpart to the access to the source code and  rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty  and the software's author,  the holder of the
% economic rights,  and the successive licensors  have only  limited
% liability. 
%
% In this respect, the user's attention is drawn to the risks associated
% with loading,  using,  modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean  that it is complicated to manipulate,  and  that  also
% therefore means  that it is reserved for developers  and  experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or 
% data to be ensured and,  more generally, to use and operate it in the 
% same conditions as regards security. 
% 
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL license and that you accept its terms.
%
%
%
% function data = createRealisticMixtures(param)
% Creates simulated mixtures of r=15 NMR spectra.
% The peaks information of the NMR spectra come from the Spectral Database for
% Organic Compounds, SDBS
% http://sdbs.riodb.aist.go.jp/sdbs/cgi-bin/cre_index.cgi
%
% Input: param structure with fields
% - m: number of observations
% - r: number of sources (max: 15)
% - Aalpha: shape parameter of the generalized gaussian distribution of A
% - Abernoulli: activation parameter of the Bernoulli distribution of A
% - width: width of the Laplacian kernel which is convolved with the peaks
% - dB: level of noise ("Inf" for noiseless data)
% - noPermute (optional): do not permute the order of the sources (default: 0)
% The number n of samples is fixed to 1024
%
%Output data structure with fields
% - A of size m x 15, absolute value of a bernoulli generalized gaussian
% matrix with activation rate A_bernoulli and shape parameter A_alpha.
% - S of size 15 x n for which is line is the spectrum of a compound.
% - names: cell-array of the names of each compound.
% - Y = A*S. A is renormalized so that norm(Y)=sqrt(m x n) (variance of
% coefficients of 1).
% - N of size m x n with standard deviation computed in order to obtain the
% specified dB level.
% - ppm: coordinates of the spectra (chemical shift)

function data = createReducedRealisticMixtures(param)

if ~isfield(param, 'width')
    param.width = 4;
end
if ~isfield(param, 'r')
    param.r = 15;
end


param.n = 2048; %spectra will are cut to 1024 size afterwards
all.spectra = zeros(18, param.n);
all.names = cell(1, 18);

 %% quinine
H1NMR400 = [3368.65   8.430    174;
   3364.14   8.418    182;
   3149.41   7.881    168;
   3140.38   7.858    182;
   2976.68   7.449    149;
   2972.17   7.437    146;
   2904.66   7.269     97;
   2901.98   7.262    124;
   2895.51   7.246     71;
   2892.94   7.239    141;
   2888.67   7.229    185;
   2886.11   7.222    124;
   2299.68   5.755     35;
   2291.99   5.735     39;
   2289.43   5.729     43;
   2282.10   5.711     64;
   2274.90   5.693     46;
   2272.34   5.686     47;
   2264.65   5.667     43;
   2193.48   5.489     86;
   2190.80   5.482     89;
   2149.66   5.379     71;
   1981.81   4.959    103;
   1964.60   4.916     99;
   1963.26   4.913     86;
   1961.91   4.910     97;
   1960.82   4.907    109;
   1950.44   4.881     95;
   1549.56   3.878   1000;
   1383.06   3.461     35;
   1382.69   3.460     35;
   1380.49   3.455     38;
   1378.30   3.449     37;
   1378.05   3.449     37;
   1225.59   3.067     85;
   1215.45   3.042     68;
   1211.79   3.033     69;
   1201.66   3.007     60;
   1055.30   2.641     53;
   1051.03   2.630     66;
   1048.58   2.624     56;
   1041.50   2.607     74;
   1039.43   2.601     72;
   1039.06   2.600     72;
   1037.72   2.597     71;
    893.07   2.235     42;
    892.82   2.235     42;
    892.46   2.234     42;
    892.09   2.233     42;
    713.75   1.786     71;
    711.55   1.781     79;
    708.50   1.773     66;
    705.81   1.767     58;
    697.02   1.745     59;
    691.89   1.732     71;
    684.69   1.714     80;
    674.80   1.689     31;
    674.44   1.688     31;
    600.95   1.504     38;
    600.46   1.503     38;
    598.02   1.497     43;
    588.38   1.473     76;
    583.01   1.459     50;
    581.18   1.455     53;
    580.81   1.454     53;
    578.86   1.449     47;
    578.37   1.448     47];
 
ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
[all.spectra(1, :), data.ppm] = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{1} = 'quinine';

%% mannitol
H1NMR400 = [   1770.87   4.432    965;
   1765.50   4.418   1000;
   1742.92   4.362    382;
   1737.18   4.347    917;
   1731.57   4.333    409;
   1663.45   4.163    772;
   1656.37   4.145    819;
   1453.25   3.637    164;
   1449.95   3.629    191;
   1447.63   3.623    172;
   1444.21   3.614    225;
   1442.50   3.610    236;
   1439.33   3.602    224;
   1436.77   3.596    215;
   1433.47   3.587    217;
   1423.34   3.562    236;
   1415.77   3.543    470;
   1408.33   3.524    357;
   1392.58   3.485    104;
   1389.04   3.476    112;
   1386.72   3.470    212;
   1383.54   3.462    212;
   1380.98   3.456    198;
   1378.05   3.449    217;
   1375.12   3.441    128;
   1372.68   3.435    110;
   1369.51   3.427    101;
   1361.69   3.408    263;
   1356.08   3.394    343;
   1350.95   3.381    275;
   1345.21   3.366    299;
   1339.23   3.352    148];

ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(2, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{2} = 'mannitol';


%% menthone (mint flavour)
H1NMR400 = [949.10   2.375    84;
    947.02   2.370     91;
    945.19   2.366    102;
    942.99   2.360     99;
    936.28   2.343    108;
    934.08   2.338    102;
    932.37   2.333    125;
    930.05   2.328    127;
    867.07   2.170     41;
    865.48   2.166     57;
    860.23   2.153     85;
    858.64   2.149     79;
    854.61   2.139     70;
    853.52   2.136    111;
    851.93   2.132     66;
    846.68   2.119     90;
    844.97   2.115     70;
    842.04   2.107     36;
    841.55   2.106     37;
    839.84   2.102     50;
    838.01   2.097     34;
    837.04   2.095     37;
    834.72   2.089     58;
    833.01   2.085     56;
    831.42   2.081     99;
    829.22   2.075     89;
    825.68   2.067    183;
    822.02   2.057    172;
    820.80   2.054    176;
    820.19   2.053    166;
    818.36   2.048    130;
    817.14   2.045    114;
    814.94   2.040     59;
    812.99   2.035     94;
    811.65   2.031    109;
    808.84   2.024    127;
    807.62   2.021    132;
    805.18   2.015     43;
    796.39   1.993    190;
    795.29   1.990    190;
    783.57   1.961    157;
    782.35   1.958    147;
    775.88   1.942     30;
    773.68   1.936     36;
    772.95   1.935     37;
    771.61   1.931     53;
    769.41   1.926     65;
    768.07   1.922     87;
    766.97   1.920     82;
    766.11   1.917     84;
    764.77   1.914     81;
    762.45   1.908     72;
    760.01   1.902     53;
    758.54   1.899     68;
    756.47   1.893     98;
    753.05   1.885     80;
    750.85   1.879     64;
    749.51   1.876     47;
    747.07   1.870     67;
    744.87   1.864     44;
    742.92   1.859     56;
    740.84   1.854     56;
    738.89   1.849     42;
    736.69   1.844     52;
    734.62   1.839     45;
    732.30   1.833     33;
    730.22   1.828     38;
    562.74   1.409     59;
    561.65   1.406     47;
    560.06   1.402     73;
    555.42   1.390     33;
    551.76   1.381    101;
    550.54   1.378    151;
    549.68   1.376    169;
    547.36   1.370    116;
    545.29   1.365     83;
    544.56   1.363     79;
    542.72   1.358    124;
    540.77   1.354     92;
    539.55   1.351     93;
    536.13   1.342     40;
    535.52   1.340     37;
    534.42   1.338     44;
    532.47   1.333     60;
    531.74   1.331     59;
    529.54   1.326     33;
    408.08   1.022    882;
    401.73   1.006    968;
    396.12   0.992     41;
    372.80   0.933     57;
    368.65   0.923    985;
    361.94   0.906    986;
    356.69   0.893     37;
    353.03   0.884     32;
    352.42   0.882     32;
    344.60   0.863    981;
    337.77   0.846   1000];

ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(3, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{3} = 'menthone';


%% caffein
H1NMR400 = [     3003.17   7.512     96;
   3002.56   7.510     95;
   1598.02   3.997    523;
   1597.41   3.996    507;
   1432.56   3.583   1000;
   1361.69   3.406    998];

ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(4, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{4} = 'caffein';


%% saccharose (sugar)
H1NMR400 = [2167.24   5.423    308;
   2163.45   5.414    314;
   1691.28   4.232    391;
   1682.50   4.210    506;
   1628.91   4.076    207;
   1620.48   4.055    344;
   1611.82   4.034    198;
   1565.80   3.918     80;
   1561.89   3.909     96;
   1559.69   3.903    166;
   1557.37   3.897     84;
   1555.79   3.893    140;
   1553.83   3.888     97;
   1551.27   3.882    135;
   1547.49   3.873    138;
   1542.85   3.861    149;
   1539.92   3.854    140;
   1537.60   3.848     99;
   1537.11   3.847     99;
   1536.50   3.845     98;
   1530.40   3.830    862;
   1527.83   3.823    854;
   1524.41   3.815    550;
   1518.07   3.799     55;
   1514.04   3.789    181;
   1504.15   3.764    267;
   1494.87   3.741    231;
   1470.34   3.680   1000;
   1430.42   3.580    235;
   1426.51   3.570    228;
   1420.41   3.555    183;
   1416.50   3.545    187;
   1397.58   3.498    200;
   1388.31   3.474    279;
   1378.78   3.450    136];

ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(5, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{5} = 'saccharose';


%% phenylalanine ("light" sugar)
H1NMR400 = [2984.01   7.467    187;
   2982.54   7.463    334;
   2980.83   7.459    165;
   2977.78   7.451    239;
   2975.71   7.446    880;
   2974.12   7.442    612;
   2969.73   7.431    406;
   2968.38   7.428    950;
   2967.53   7.426    672;
   2966.19   7.422    207;
   2962.77   7.414     50;
   2960.69   7.409    214;
   2959.23   7.405    475;
   2957.64   7.401    365;
   2955.81   7.396     67;
   2954.47   7.393    143;
   2952.03   7.387    514;
   2948.85   7.379    148;
   2946.04   7.372    141;
   2944.46   7.368    231;
   2942.50   7.363    823;
   2940.92   7.359   1000;
   2938.60   7.353    251;
   2935.91   7.347    348;
   2934.33   7.343    640;
   1606.08   4.019    471;
   1600.83   4.006    525;
   1598.14   3.999    536;
   1595.34   3.992     52;
   1592.90   3.986    505;
   1334.23   3.339    221;
   1328.98   3.326    216;
   1319.58   3.302    328;
   1314.33   3.289    309;
   1260.13   3.154    365;
   1252.08   3.133    351;
   1245.48   3.117    252;
   1237.43   3.097    236];

ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(6, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{6} = 'phenylalanine';




%% lactose (milk)
H1NMR400=[2534.55   6.342    670;
   2530.15   6.331    662;
   2040.77   5.107    510;
   2036.99   5.097    510;
   1961.30   4.908    407;
   1957.15   4.898    618;
   1953.13   4.888    420;
   1915.41   4.793    354;
   1912.23   4.785    361;
   1880.98   4.707     45;
   1872.07   4.685    345;
   1865.97   4.670    565;
   1861.82   4.659    606;
   1856.57   4.646    292;
   1849.00   4.627     43;
   1846.44   4.621     55;
   1823.97   4.564     41;
   1823.73   4.564     41;
   1812.87   4.537    460;
   1808.35   4.525    468;
   1801.15   4.507     79;
   1795.78   4.494    262;
   1789.67   4.479    562;
   1785.03   4.467   1000;
   1680.42   4.205     41;
   1679.81   4.204     42;
   1679.57   4.203     42;
   1679.20   4.202     43;
   1672.24   4.185    481;
   1664.79   4.166    480;
   1659.79   4.154     41;
   1657.59   4.148     41;
   1488.53   3.725    175;
   1484.99   3.716    368;
   1481.69   3.708    244;
   1478.64   3.700    206;
   1475.10   3.691    390;
   1471.92   3.684    282;
   1468.02   3.674     66;
   1466.55   3.670     69;
   1457.52   3.647    593;
   1452.03   3.634    790;
   1448.36   3.625    840;
   1436.04   3.594     71;
   1434.33   3.589     64;
   1429.81   3.578    284;
   1425.17   3.567    142;
   1420.65   3.555    668;
   1414.92   3.541    361;
   1411.01   3.531    560;
   1404.91   3.516    362;
   1397.71   3.498    387;
   1392.46   3.485    326;
   1386.84   3.471    216;
   1382.20   3.459    721;
   1376.46   3.445    457;
   1370.12   3.429    244;
   1363.65   3.413    115;
   1360.84   3.406    101;
   1359.25   3.402     95;
   1358.64   3.400     93;
   1357.30   3.397     90;
   1355.96   3.393     85;
   1355.22   3.392     84;
   1354.86   3.391     82;
   1354.13   3.389     80;
   1353.52   3.387     79;
   1351.93   3.383     76;
   1350.46   3.380     73;
   1349.37   3.377     72;
   1348.88   3.376     71;
   1348.27   3.374     69;
   1348.02   3.374     69;
   1347.66   3.373     69;
   1347.41   3.372     69;
   1346.80   3.370     68;
   1346.07   3.369     68;
   1345.70   3.368     67;
   1345.34   3.367     67;
   1344.97   3.366     67;
   1344.36   3.364     66;
   1343.99   3.363     66;
   1343.51   3.362     65;
   1342.90   3.361     64;
   1340.09   3.354    109;
   1335.94   3.343     90;
   1330.44   3.330    351;
   1323.85   3.313    816;
   1312.62   3.285    495;
   1303.71   3.263    564;
   1294.31   3.239    331;
   1276.12   3.194    148;
   1270.51   3.180    215;
   1267.09   3.171    257;
   1263.55   3.162    198;
   1257.81   3.148    135];

ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(7, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{7} = 'lactose';


%% ascorbic acid (C vitamin)
H1NMR400 = [3320.62   8.309     81;
   1940.92   4.857     52;
   1912.54   4.786     32;
   1911.93   4.785     32;
   1911.16   4.783     31;
   1910.71   4.781     31;
   1887.51   4.723    988;
   1885.99   4.720   1000;
   1500.09   3.754    133;
   1498.72   3.751    134;
   1492.31   3.735    321
   1490.94   3.731    222;
   1485.90   3.719    167;
   1484.53   3.715    169;
   1389.47   3.477     66;
   1383.06   3.461    103;
   1378.94   3.451    561;
   1377.11   3.446    533;
   1372.68   3.435    432;
   1369.17   3.426    429;
   1366.58   3.420    118;
   1358.64   3.400     77];

ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(8, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{8} = 'ascorbic acid';


%% citric acid (lemon)
H1NMR400 = [1112.06   2.783    570;
   1096.68   2.745   1000;
   1069.58   2.677    994;
   1054.20   2.638    559];

ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(9, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{9} = 'citric acid';


%% sodium p-hydroxybenzoate (in Orangina)
H1NMR400 = [3120.24   7.808    919;
   3111.82   7.787   1000;
   2694.75   6.743    991;
   2686.14   6.722   1000];

ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(10, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{10} = 'sodium p-hydroxybenzoate';


%% alcohol in fruits
H1NMR400 = [4963.26  12.420    116;
   1712.79   4.286    727;
   1708.21   4.275    822;
   1705.10   4.267    892;
   1700.33   4.255    770;
   1058.15   2.648    658;
   1053.39   2.636    684;
   1042.40   2.609    995;
   1037.64   2.597    970;
    991.48   2.481   1000;
    983.60   2.462    970;
    975.91   2.442    660;
    968.03   2.423    661];

ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(11, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{11} = 'fruit alcohol';


%% folic acid
H1NMR400=[3467.02   8.676   1000;
   3271.50   8.186    313;
   3263.81   8.167    315;
   3073.20   7.690    717;
   3064.44   7.668    761;
   2787.51   6.975    333;
   2669.56   6.680    726;
   2660.79   6.658    722;
   1806.91   4.522    519;
   1801.36   4.508    514;
   1749.00   4.377    138;
   1746.01   4.369    164;
   1744.52   4.366    184;
   1741.53   4.358    174;
   1739.82   4.354    168;
   1739.39   4.353    168;
   1736.82   4.346    147;
    944.91   2.365    301;
    937.65   2.347    655;
    930.17   2.328    374;
    832.73   2.084    139;
    825.03   2.065    165;
    773.96   1.937    143];

ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(12, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{12} = 'folic acid';


%% taurine
H1NMR400 = [1525.76   3.818   1000;
   1446.04   3.619     80;
   1439.33   3.602    181;
   1432.74   3.585    102;
   1339.48   3.352    106;
   1332.89   3.336    195;
   1326.17   3.319     81];

ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(13, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{13} = 'taurine';

   
%% cholesterol 
H1NMR400 = [2140.26   5.356     88;
   2138.43   5.351     67;
   2137.08   5.348     65;
   2135.01   5.343     81;
   1413.70   3.538     31;
   1412.23   3.534     45;
   1408.08   3.524     54;
   1407.10   3.521     55;
   1402.83   3.511     45;
   1401.25   3.507     35;
   1396.00   3.494     30;
    907.10   2.270     79;
    797.49   1.996     99;
    740.72   1.854    154;
    731.32   1.830    157;
    596.68   1.494    183;
    530.88   1.329    124;
    439.45   1.100    168;
    402.22   1.007   1000;
    368.77   0.923    498;
    362.18   0.907    465;
    349.73   0.876    764;
    347.90   0.871    757;
    343.02   0.859    717;
    341.31   0.855    699;
    276.37   0.692     42;
    270.87   0.678    927];
    
ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(14, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{14} = 'cholesterol';

    
%% adenosine
H1NMR400 = [3349.61   8.382    942;
   3269.90   8.182     38;
   3269.53   8.181     39;
   3264.28   8.168   1000;
   2960.69   7.409    448;
   2363.77   5.915    431;
   2357.54   5.900    439;
   2202.03   5.510    388;
   2195.80   5.495    475;
   2191.16   5.483    194;
   2186.52   5.472    134;
   2098.02   5.250    333;
   2093.51   5.239    330;
   1863.53   4.663    103;
   1857.54   4.648    221;
   1852.29   4.635    222;
   1846.19   4.620    101;
   1673.22   4.187    111;
   1668.70   4.176    215
   1665.77   4.169    218;
   1661.13   4.157    109;
   1600.83   4.006    119;
   1597.66   3.998    304;
   1594.48   3.990    296;
   1591.19   3.982    118;
   1487.92   3.724     67;
   1484.25   3.714    113;
   1480.35   3.705     80;
   1475.71   3.693    120;
   1474.61   3.690    110;
   1472.05   3.684    178;
   1468.63   3.675    101;
   1441.16   3.607     96;
   1437.62   3.598    114;
   1434.33   3.589    115;
   1433.96   3.589    115;
   1430.18   3.579    128;
   1425.66   3.568     79;
   1422.00   3.559     77;
   1418.46   3.550     63];
   
ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(15, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{15} = 'adenosine';

   
%% myo-inositol (inositol: in red bull)
H1NMR400 = [1627.44   4.073    326;
   1624.63   4.066    661;
   1621.70   4.058    401;
   1459.96   3.654    465;
   1450.07   3.629   1000;
   1444.46   3.615     43;
   1440.67   3.605    919;
   1428.10   3.574     22;
   1427.12   3.571     23;
   1425.78   3.568     25;
   1425.17   3.567     25;
   1420.41   3.555    926;
   1417.60   3.548    938;
   1410.40   3.530    507;
   1407.59   3.523    554;
   1320.80   3.305    414;
   1311.52   3.282    641;
   1302.25   3.259    306];
   
ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(16, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{16} = 'myo-inositol';

   
%% oleic acid (unsaturated fac)
H1NMR400 = [2142.46   5.361    102;
   2141.36   5.359     81;
   2138.79   5.352    131;
   2136.72   5.347    237;
   2135.74   5.345    243;
   2133.67   5.339    138;
   2133.18   5.338    138;
   2130.98   5.333     81;
   2130.00   5.330    109;
    944.95   2.365    171;
    937.50   2.346    306;
    929.81   2.327    202;
    806.76   2.019    261;
    800.90   2.005    246;
    794.92   1.990    111;
    659.42   1.650    101;
    652.34   1.633    150;
    645.02   1.614    127;
    639.04   1.599     61;
    637.57   1.596     62;
    523.93   1.311    812;
    515.01   1.289    540;
    506.71   1.268   1000;
    358.52   0.898    250;
    351.81   0.881    774;
    344.73   0.863    303];
    
ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(17, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{17} = 'oleic acid';


%% glycerol (coca)
H1NMR400 = [1792.36   4.485    411;
   1390.63   3.480     72;
   1384.89   3.466    238;
   1379.64   3.453    453;
   1374.39   3.439    346;
   1368.90   3.426    203;
   1360.72   3.405    624;
   1355.59   3.392    316;
   1349.85   3.378   1000;
   1344.73   3.365    712;
   1334.72   3.340     34;
   1329.35   3.327    993;
   1323.49   3.312    744;
   1318.48   3.300    574;
   1312.62   3.285    435];

ppm = H1NMR400(:, 2)';
peak = H1NMR400(:, 3)';
all.spectra(18, :) = Get_NMR_Spectrum(peak, ppm, [0, 10], param.n, param.width);
all.names{18} = 'glycerol';


%% Create source and mixture matrices
% sources
indices = [3, 4, 6, 7, 8, 9, 13, 11, 12, 2, 14, 5, 16, 17, 18];
data.S = all.spectra(indices, :);
data.nS = all.spectra(indices, :);
data.names = all.names(indices);


% reduce the size and the number
param.n = 1024;
range = 112 : 1135;
perm = 1 : 15;
if isfield(param, 'noPermute')
    if param.noPermute == 0
        perm = randperm(15);
    end
end
perm = perm(1 : param.r);
data.names = data.names(perm);
data.S = data.S(perm, range);
data.nS = data.nS(perm, range);
data.ppm = data.ppm(range);


% mixtures
data.A = abs(gengau2(param.Aalpha, param.m, param.r)...
    .* (rand(param.m, param.r) < param.Abernoulli));


% normalize to have std(Y_ij)=1
data.Y = data.A * data.S;
coeff = sqrt(param.m * param.n) / norm(al(data.Y), 'fro');
data.A = data.A * coeff;
data.Y = data.Y * coeff;


%% additive and multiplicative noise
% compute additive noise level
if isfinite(param.dB)
    noise = 10.^(-param.dB / 20)...
        / sqrt(numel(data.Y)) * norm(data.Y(:), 'fro');
else
    noise = 0;
end
data.N = randn(size(data.Y)) * noise;

% multiplicative noise if specified
if isfield(param, 'multStd')
    data.Nadd = data.N;
    data.Nmult = param.multStd * randn(size(data.Y)) .* data.Y;
    data.N = data.Nadd + data.Nmult;
end

end




