%Name:		Chris Shoemaker
%Course:	EER-280 - Digital Watermarking
%Project: 	Determination of Periodicity of MATLAB's rand() function
%           If a pattern is repeated, correlation should jump to 1

clear all;

% read in key for PN generator
file_name='_key.bmp';
key=double(imread(file_name))./256;

% reset MATLAB's PN generator to state "key"
rand('state',key);

% generate the baseline PN sequence
pn_base_sequence=round(2*(rand(8,8)-0.5));

% go up to 120,000 patterns
for kk=1:120000
    
    % display progress to the console
    if (mod(kk,10000) == 0)
        kk,
    end    
    
    % generate a new sequence and calulate correlation to baseline
    pn_sequence=round(2*(rand(8,8)-0.5));
    correlation(kk)=corr2(pn_sequence,pn_base_sequence);
end

% plot the correlation for each sequence
figure(1)
set(1,'color','white')
plot(correlation(1:kk))
title('Correlation Between PN Sequences Generated using MATLAB rand() vs Baseline Sequence')
xlabel('PN Sequence')
ylabel('Correlation')
