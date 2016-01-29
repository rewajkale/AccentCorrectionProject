function out=psola(in,m,alpha,beta)
% in input signal
% m pitch marks (from PitchMarker.m function)
% alpha time stretching factor
% beta pitch shifting factor
P = diff(m); %compute pitch periods
if m(1)<=P(1), %remove first pitch mark
m=m(2:length(m));
P=P(2:length(P));
end
if m(length(m))+P(length(P))>length(in) %remove last pitch mark
m=m(1:length(m)-1);
else
P=[P P(length(P))];
end
Lout=ceil(length(in)*alpha);
out=zeros(1,Lout); %output signal
tk = P(1)+1; %output pitch mark
while round(tk)<Lout
[minimum i] = min( abs(alpha*m - tk) ); %find analysis segment
pit=P(i);
st=m(i)-pit;
en=m(i)+pit;
if st<0
    gr=in(st-st+1:en-st+1) .* hanning(2*pit+1);
else
gr = in(st:en) .* hanning(2*pit+1);
end
iniGr=round(tk)-pit;
endGr=round(tk)+pit;
if endGr>Lout, break; 
end
if iniGr<0
    out(iniGr-iniGr+1:endGr-iniGr+1) = out(iniGr-iniGr+1:endGr-iniGr+1)+gr';
else
out(iniGr:endGr) = out(iniGr:endGr)+gr'; %overlap new segment
end
tk=tk+pit/beta;
end
end