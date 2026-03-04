function [formants, heights, index,BW] = formant_bw_id( spectra, freq,numformants);
%Author: Brad Story
%Modified: 11.03.2025
% peak detector looks for any positive maxima and finds bandwidths
% Used in spectrum formant tracking, since spectra is
% smoothly varying. Assumes a cascade generated spectra 

%INPUTS: 
% spectra = 20log10(abs(h)) if using with VokalTrakM.m
% freq = frequency vector that coincides with "spectra"
% numformants = number of formants the algorithm should return

%OUTPUTS:
% formants = frequencies of the peaks in spectra
% heights = amplitude of the peaks in spectra
% index = the freq vector index of the peaks in spectra
% BW = estimated bandwidths of each peak in spectra

%Example call (if coming after using VokalTraktM.m):
%[fmts, hts, index,BW] = formant_bw_id( 20*log10(abs(h)),f,4);

siglen = length(spectra) ;
; 
formants = [] ;
heights = [] ;
index = [] ;
df = freq(2) - freq(1);
for i = 2:siglen-1
   if( (spectra(i-1) < spectra(i)) & ...
         (spectra(i+1) < spectra(i) ) )
      fprm = (spectra(i+1) - spectra(i-1))/(2*df);
      fprm2 = (spectra(i+1) - 2*spectra(i)+spectra(i-1))/(df^2);
      corr = -fprm/fprm2;
      formants = 	[formants  freq(i)+corr ] ;
      heights  =  [heights   spectra(i)] ; 
      index = [index i];
   end ; 
end ; 

if (length(formants) > numformants)
   formants = formants(1:numformants);
   heights = heights(1:numformants);
end;




%========Bandwidths===============

f = freq;
Df = f(2)-f(1); 


for n=1:length(formants)
    Stmp = spectra-heights(n)+3;
  
    i = index(n);
    while( Stmp(i) > 0)
        i=i-1;
        if(i<=1)
            i=1;
            break;
        end
    end
    indexL = i;
    
    df = -(Df*Stmp(i))/(Stmp(i+1) - Stmp(i));
    fL(n) = f(i)+df;
       
    i = index(n);
    while( Stmp(i) > 0)
        i=i+1;
        if(i>=length(Stmp))
            i=length(Stmp)-1;
            break;
        end
    end
    indexR = i;
    
    df = -(Df*Stmp(i))/(Stmp(i+1) - Stmp(i));
    fR(n) = f(i)+df;
   
    fR(n);
    fL(n);
   BW(n) = fR(n)-fL(n); 
   
    
end


    