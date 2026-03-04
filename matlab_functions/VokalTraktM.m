function [h,z,f,D,h_alt,Zrad,z_alt] = VokalTraktM(ARvect,LTvect,cutoff,maxF)
%Cascade (transmission line) calculation of vocal tract frequency response
%and input impedance
%Author: B. Story (origination May 1995, revised 11.02.25)
%Based on Sondhi and Schroeter 1987 and Story et al. 2000
%
%INPUTS
%ARvect = vector of cross-sectional areas representing a vocal tract area
%function. Element "1" of the vector is located just above the glottis.
%This vector can have any number of elements but must be greater than one
%element
%LTvect = vector of section/tublet lengths. Must be the same number of
%elements as ARvect. Lengths can be the same or different
%cutoff = this is number of sections for which a yielding wall will be
%imposed. If hard walled tubes are desired set cutoff=0; if yielding walls
%are desired set cutoff = number of elements in ARvect + 1
%maxF = the maximum frequency for which to calculate the frequency response
%and input impedance
%OUTPUTS:
%h = complex vector containing the frequency response. 
%z = complex vector containing the input impedance
%f = frequency vector
%h_alt = complex vector containing the frequency response calculated
%without radiation impedance
%Zrad = radiation impedance
%z_alt = complex vector containing the input impedance calculated w/o
%radiation impedance
%
% EXAMPLE CALL:  [h,z,f,D,h_alt] = VokalTraktM(ones(44,1),0.396825*ones(44,1),45,5000);
% or [h,z,f,D,h_alt] = VokalTraktM(arvect,dxvect,50,5000); 
% where arvect = 44-element vector of cross-sectional areas
% and dxvect = 44-element vector of length for each section
%EXAMPLE PLOT:
%to plot the magnitude of the frequency response: >>plot(f,20*log10(abs(h)));
%to plot the magnitude of the input impedance: >>plot(f,20*log10(abs(z)));

NSect 		= length(ARvect);
NPoints		= round(maxF/1);				% No. of points in Freq Spect
%df 			= round(maxF/NPoints);				% delta f - freq interval
df 			= (maxF/NPoints);

r = 408;
Ft = 250;
%Use Ft = 340 for female vocal tract area functions
%disp('Now running the yielding wall for females!!');
%Ft = 340;
Fw = 15;

q 			= 4.0;
csnd 		= 35000.0;				% speed of sound
rho 		= 0.00114;				% density of air
rhoc		= rho*csnd;
FreqVect 	= [df:df:df*NPoints]';

FreqVectLT	= length(FreqVect);

om 			= 2*pi*FreqVect;
alpha 		= zeros(FreqVectLT,1) + j*(q*om);
alpha 		= sqrt(alpha);
temp1 		= r*ones(FreqVectLT,1) + j*(om);
temp2 		= zeros(FreqVectLT,1) + j*(om);
den 		= temp1.*temp2;
temp1 		= ((2*pi*Fw)^2)*ones(FreqVectLT,1) + j*zeros(FreqVectLT,1);
den 		= den + temp1;
num 		= zeros(FreqVectLT,1) + j*om*(2*pi*Ft)^2;
beta 		= (num ./ den) + alpha;


%disp('Number of Vocal Tract Sections: ');
%NSect

A 			= ones(NPoints,1);
B 			= zeros(NPoints,1);
C 			= zeros(NPoints,1);
D 			= ones(NPoints,1);

for k=1:NSect
  
  
  %----------This is for adding a hard-walled tube to the end of the tract --
  
  if(k >=  cutoff)

    r = 0;
    beta = 0;
    alpha = 0;
    
  end;


  temp1 		= r*ones(FreqVectLT,1) + j*(om);
  temp2 		= beta + (zeros(FreqVectLT,1) + j*om);
  gamma 		= sqrt(temp1 ./ temp2);
  sigma 		= gamma .* temp2;


% -------------------------------------------------------------------  

  
	nA 		= cosh( (LTvect(k)/csnd) * sigma);
	temp1 	= -rho*csnd/ARvect(k) * gamma;
	nB 		= temp1 .* sinh( (LTvect(k)/csnd) * sigma);
	temp1 	= ones(length(gamma),1) ./ gamma;
	temp2 	= -ARvect(k)/(rho*csnd) * temp1;
	nC 		= temp2 .* sinh((LTvect(k)/csnd) * sigma);
	nD 		= nA;
	
  
        
        M1 		= [A B; C D];
		M2 		= [nA nB; nC nD];
		

        A      = [M2(1:NPoints,1).*M1(1:NPoints,1) + M2(1:NPoints,2).*M1(NPoints+1:end,1)];
        B      = [M2(1:NPoints,1).*M1(1:NPoints,2) + M2(1:NPoints,2).*M1(NPoints+1:end,2)];
        C      = [M2(NPoints+1:end,1).*M1(1:NPoints,1) + M2(NPoints+1:end,2).*M1(NPoints+1:end,1)];
        D      = [M2(NPoints+1:end,1).*M1(1:NPoints,2) + M2(NPoints+1:end,2).*M1(NPoints+1:end,2)];
end



  R     = 128. * rho*csnd/(9*pi^2 * ARvect(NSect));    
  L		= 8.*rho*csnd/((3*pi*csnd)*sqrt(ARvect(NSect)*pi)); 	
    

%Calculate radiation impedance  
temp1 	= zeros(NPoints,1) + j*R*L*om;
temp2 	= R*ones(NPoints,1) + j*L*om;
Zrad	= temp1 ./ temp2;

%Calculate input impedance z and freq response h
temp1	= (Zrad .* D) - B;
temp2 	= A - (C .* Zrad);
z		= temp1 ./ temp2;
h		= Zrad ./ temp2;
f = FreqVect;


%Calculate h and z without radiation
h_alt = (-C.*B)./A +D;
z_alt = -B./A;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


