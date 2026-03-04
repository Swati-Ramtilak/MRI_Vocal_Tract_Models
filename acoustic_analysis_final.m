%------------ Load Cross sectional areas of Conventional GRE and Hybrid Models-----------------
clear;clc;
try1 = 'enter path to gre vocal tract area function';   
csvdata_gre = readtable(try1);
dxvect_gre = csvdata_gre.Distance_from_L0_cm.'; 
save('dxvect_gre.mat', 'dxvect_gre');
aa_gre     = csvdata_gre.Cross_Sectional_Area_cm2.'; 
save('aa_gre.mat', 'aa_gre');

try2 = 'enter path to gre+teeth vocal tract area function';   
csvdata_teeth = readtable(try2);
dxvect_zte = csvdata_teeth.Distance_from_L0_cm.'; 
save('dxvect_teeth.mat', 'dxvect_teeth');
aa_teeth     = csvdata_teeth.Cross_Sectional_Area_cm2.'; 
save('aa_teeth.mat', 'aa_teeth');

%------------ Calculate freq response of aa_gre-----------------
load('aa_gre.mat')
load('dxvect_gre.mat')
dx = dxvect_gre(end)/40; %find the uniform section length
xnew = [dx:dx:dxvect_gre(end)]; %generate new cummulative length vector
aa_gre_new = interp1(dxvect_gre,aa_gre,xnew,'pchip'); %resample area vector based on constant section length
aa_gre_new = aa_gre_new(:); %make sure this is a column vector
%Generate section length vector - typically all entries are the same
%(i.e.,dx)
dxvect = dx*ones(44,1); %sets the section length vector based on the dx calculated above
%                        
cutoff = 50; %this assures that the entire vocal tract will have energy losses
maxF = 5000; %maximum frequency in calculation - typically 5000 Hz, but can be increased or decreased
[h1,z,f1,D,h_alt,Zrad] = VokalTraktM(aa_gre_new,dxvect,cutoff,maxF);
[fmts1, hts, index,BW] = formant_bw_id( 20*log10(abs(h1)),f1,4);

%------------ Calculate freq response of aa_zte-----------------
load('aa_zte.mat')
load('dxvect_zte.mat')
dx = dxvect_gre(end)/40; %find the uniform section length
xnew = [dx:dx:dxvect_gre(end)]; %generate new cummulative length vector
aa_zte_new = interp1(dxvect_zte,aa_teeth,xnew,'pchip'); %resample area vector based on constant section length
aa_zte_new = aa_zte_new(:); %make sure this is a column vector
%Generate section length vector - typically all entries are the same
%(i.e.,dx)
dxvect = dx*ones(44,1); %sets the section length vector based on the dx calcul
cutoff = 50; %this assures that the entire vocal tract will have energy losses
maxF = 5000; %maximum frequency in calculation - typically 5000 Hz, but can be increased or decreased
[h2,z,f2,D,h_alt,Zrad] = VokalTraktM(aa_zte_new,dxvect,cutoff,maxF);
[fmts2, hts, index,BW] = formant_bw_id( 20*log10(abs(h2)),f2,4);

%-----Plot results from gre and zte-----------------
figure(1)
clf; hold on
p1 = plot(f1(1:5:end),20*log10(abs(h1(1:5:end))),'-b','LineWidth',1);
p2 = plot(f2(1:5:end),20*log10(abs(h2(1:5:end))),'-r','LineWidth',1);
for n=1:length(fmts1)
    plot([fmts1(n) fmts1(n)],[-50 100],'-b','LineWidth',0.5);
end
for n=1:length(fmts2)
    plot([fmts2(n) fmts2(n)],[-50 100],'-r','LineWidth',0.5);
end
axis([0 maxF -20 50]); %change as needed
set(gca,'PlotBoxAspectRatio',[1 .5 1],'FontSize',12);
xlabel('Frequency (Hz)')
ylabel('Rel. ampl. (dB)')
grid on
hL = legend([p1 p2],{'GRE-MRI','GRE-MRI + ZTE Teeth'});