%Name:		Chris Shoemaker
%Course:	EER-280 - Digital Watermarking
%Project: 	Tests that using rand X times generates the same values as rand(1,X)

clear all;

% read in key for PN generator
file_name='_key.bmp';
key=double(imread(file_name))./256;

% reset MATLAB's PN generator to state "key"
rand('state',key);

% execute rand 20 times
for i=1:20
    test1(i)=rand;
end

% reset MATLAB's PN generator to state "key"
rand('state',key);

% generate 20 element random vector
test2=rand(1,20);

% check to make sure they're equal
difference=sum(test1-test2),