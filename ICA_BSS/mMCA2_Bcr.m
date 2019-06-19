function [part,options]=mMCA2_Bcr(X,dict,pars1,pars2,pars3,itermax,gamma,comptv,expdecrease,stop,mask,sigma,display,numic)

% Dictionary metadata.
numberofdicts = LengthList(dict);
options.nbcomp= sprintf('Number of morphological components: %d',numberofdicts);
options.dict  = ['Transforms: [ ' dict ']'];
str = 'Parameter 1 of transforms: [ ';
for nb=1:numberofdicts, str = [str num2str(NthList(pars1,nb)) ' : ']; end
options.pars1 = [str ']'];
str = 'Parameter 2 of transforms: [ ';
for nb=1:numberofdicts, str = [str num2str(NthList(pars2,nb)) ' : ']; end
options.pars2 = [str ']'];	
str = 'Parameter 3 of transforms: [ ';
for nb=1:numberofdicts, str = [str num2str(NthList(pars3,nb)) ' : ']; end
options.pars3 = [str ']'];
part	      = zeros(numic,65536,numberofdicts);	

A = rand(numic,numic);
part = X;
imgpad = X;
% Start the modified Block Relaxation Algorithm.
for iter=0:itermax-1
	%for i=1:J
	  % Calculate the residual image.
	    residual=imgpad-sum(A*part);
	  % Cycle over dictionaries.
	   for nb=1:numberofdicts
	   % Update Parta assuming other parts fixed.
	   % Solve for Parta the marginal penalized minimization problem (Hard thesholding, l_1 -> Soft).
	     NAME   = NthList(dict,nb);
	     PAR1   = NthList(pars1,nb);
	     PAR2   = NthList(pars2,nb);
	     PAR3   = NthList(pars3,nb);
         
	     Ra=part(:,:,nb)+residual;
	     coeffa = FastLA2(Ra,NAME,PAR1,PAR2,PAR3); % get coefficients c = \Phi^T * x
	     aaa;
         coeffa = eval([thdtype 'ThreshStruct(coeffa,delta,NAME);']); %soft/hard thresholdng method
	     
         part(:,:,nb)   = FastLS2(coeffa,NAME,PAR1,PAR2,PAR3);
	     if (nb == comptv) & gamma~=0, part(:,:,nb)  = TVCorrection(part(:,:,nb),gamma); end
	   end
	%end
	
	% Update the regularization parameter delta.
	    if expdecrease	delta=delta*lambda; % Exponential decrease.	
	    else		delta=delta-lambda; % Linear decrease
        end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = HardThreshStruct(C,lambda,nameofdict)
global E

nbdicts = length(C);

for nb=1:nbdicts
 coeffs = C{nb};
 scaleindex = length(coeffs);
 if strcmp(nameofdict,'CURVWRAP')
   for j = 2:scaleindex
     for w = 1:length(coeffs(j).coeff)
       coeffs(j).coeff{w} = coeffs(j).coeff{w}.* (abs(coeffs(j).coeff{w}) > lambda*E{j}{w});
     end
   end
 else
   for j = 2:scaleindex
  	  coeffs(j).coeff = coeffs(j).coeff .* (abs(coeffs(j).coeff) > lambda);
   end
 end
 C{nb} = coeffs;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = SoftThreshStruct(C,lambda,nameofdict)
global E

nbdicts = length(C);

for nb=1:nbdicts
 coeffs = C{nb};
 scaleindex = length(coeffs);
 if strcmp(nameofdict,'CURVWRAP')
   for j = 2:scaleindex
     for w = 1:length(coeffs(j).coeff)
       coeffs(j).coeff{w} = sign(coeffs(j).coeff{w}) .* max(abs(coeffs(j).coeff{w}) - lambda*E{j}{w},0);
     end
   end
 else
   for j = 2:scaleindex
  	  coeffs(j).coeff = sign(coeffs(j).coeff) .* max(abs(coeffs(j).coeff) - lambda,0);
   end
 end
 C{nb} = coeffs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delta = StartingPoint(C,dict)
global E

nbdicts = length(C);

for nb=1:nbdicts
 tmp = [];
 coeffs = C{nb};
 scaleindex = length(coeffs);
 
 % If it is curvelet basis by the wrapping algorithm, then compute the L2 norm of the basis elements
 if strcmp(NthList(dict,nb),'CURVWRAP')
   computeL2norm(coeffs);
   for j = 2:scaleindex
     for w = 1:length(coeffs(j).coeff)
       wedge = coeffs(j).coeff{w}/E{j}{w};
       tmp = [tmp;wedge(:)];
     end
   end
 else
   for j = 2:scaleindex
  	  tmp = [tmp;coeffs(j).coeff(:)];
   end
 end
 buf(nb)=max(abs(tmp(:)));
end

%
buf=flipud(sort(buf(:),1))';
if nbdicts>1 delta=buf(2);
else	     delta=buf(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = TVCorrection(x,gamma)
% Total variation implemented using the approximate (exact in 1D) equivalence between the TV norm and the l_1 norm of the Haar (heaviside) coefficients.

[n,J] = quadlength(x);

qmf = MakeONFilter('Haar');

[ll,wc,L] = mrdwt(x,qmf,1);

wc = SoftThresh(wc,gamma);

y = mirdwt(ll,wc,qmf,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function computeL2norm(coeffs)
% Compute norm of curvelets (exact)
global E

%F = ones(size(coeffs(end).coeff{1}));
F = ones(coeffs(1).coeff{2});
X = fftshift(ifft2(F)) * sqrt(prod(size(F))); 	% Appropriately normalized Dirac
C = fdct_wrapping(X,1,length(coeffs)); 		% Get the curvelets

E = cell(size(C));
for j=1:length(C)
  E{j} = cell(size(C{j}));
  for w=1:length(C{j})
    A = C{j}{w};
    E{j}{w} = sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
  end
end
