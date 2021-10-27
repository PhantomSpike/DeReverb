function y = cosrampenv(len,ramplen,fs)
%y = cosrampenv(len,ramplen,fs)
%creates an 1xN envelope y with an onset cosine ramp and an offset cosine ramp
%and a flat middle. If there is a third optional argument fs (the sampling
%rate), then the length of the envelope len, and the length of the ramps
%ramplen, are taken to be in seconds. Otherwise they are taken to be in
%samples.
%
%Nicol Harper April 25th 2016

if nargin<3
    %units are in samples
    slen = len;
    sramlen = ramplen;
else
    %units are in seconds
    slen = round(fs.*len);
    sramlen = round(fs.*ramplen);
end


x = 0:(sramlen-1);
frontramp = 0.5-0.5.*cos(pi.*x./(sramlen-1));
backramp = fliplr(frontramp);
y = [frontramp ones(1,slen-2.*sramlen) backramp];
