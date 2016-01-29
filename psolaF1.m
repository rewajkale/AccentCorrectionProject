function out=psolaF1(in,m,alpha,beta,gamma)

% gamma newFormantFreq/oldFormantFreq

P = diff(m);%compute pitch periods
pitch=max(m);
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
[minimum i]=min(abs(alpha*m-tk) ); % find analysis segment
pit=P(i);
pitStr=floor(pit/gamma);
k2=m(i)-pit;
k3=m(i)+pit;
if k2<0
    %size(in(k2-k2+1:k3-k2+1))
    %size(hanning(2*pit+1))
    gr=in(k2-k2+1:k3-k2+1).*hanning(2*pit+1);
else
    %size(in(k2:k3))
    %size(hanning(2*pit+1))
   gr=in(k2:k3).*hanning(2*pit+1); 
end
%gr=in(m(i)-pit:m(i)+pit).*hanning(2*pit+1);
gr=interp1(-pit:1:pit,gr,-pitStr*gamma:gamma:pit);% stretch segm.
iniGr=round(tk)-pitStr;endGr=round(tk)+pitStr;
if endGr>Lout, break; 
end
if iniGr<0
    out(iniGr-iniGr+1:endGr-iniGr+1)=out(iniGr-iniGr+1:endGr-iniGr+1)+gr;
else
out(iniGr:endGr)=out(iniGr:endGr)+gr; % overlap new segment
end
tk=tk+pit/beta;
end