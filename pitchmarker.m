function [ pitch ] = pitchmarker(blk_section)
blk_size=400;
mark=[1:length(blk_section)]*0;
last_pos=1;
place=1;
blk_size=300;
i=1;
while last_pos+floor(blk_size*1.7) < length(blk_section)
% grabs the next block to examine
temp=blk_section(last_pos+50:last_pos+floor(blk_size*1.7));
% finds the high point in the block
[mag,place]=max(temp);
% check for a signal in the current block
if mag < 0.01
place=length(temp);
mode = 0;
mark(place+last_pos+50)=1;
pitch(i)=place+last_pos+50;
else
mode = 1;
end
% check for pitch mark before current pitchmark
while mode == 1
% find the largest point in block from start to current pitch mark
[mag2,place2]=max(temp(1:place-50));
% check if high mark has great enough magnitude to be a pitch mark
if mag2 > 0.90*mag
mag=mag2;
place=place2;
else
mode = 0;
mark(place+last_pos+50)=1;
pitch(i)=place+last_pos+50;
end
end
% next block to look at is 50 samples after current block
blk_size=place+50;
% makes sure next blk_size is of large enough size
if blk_size < 150
blk_size=150;
end
last_pos=place+last_pos+50;
i=i+1;
end
figure(1)
hold on
plot(mark)
plot(blk_section,'r')
title('Pitch Marking Waveform')
xlabel('time')
ylabel('my.wav(t)')