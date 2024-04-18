function [sig_src,sig_rec,sig_out] = Mini_Alignment(sig_src, sig_rec)
%% Normelize RMS For sig_rec
% sig_rec = sig_rec/rms(sig_rec)*rms(sig_src);

%% Do Xcorr
if length(sig_rec) < length(sig_src)
%     [C,lags] = xcorr(abs(sig_rec),abs(sig_src));    %sig_rec is ref, sig_src is being moved upon sig_rec.
    [C,lags] = xcorr(abs(sig_src),abs(sig_rec));    %sig_sec is ref, sig_rec is being moved upon sig_src.
else
%     [C,lags] = xcorr(abs(sig_src),abs(sig_rec));    %sig_sec is ref, sig_rec is being moved upon sig_src.
    [C,lags] = xcorr(abs(sig_rec),abs(sig_src));    %sig_rec is ref, sig_src is being moved upon sig_rec.
end
Cm = abs(C(lags > 0 & lags < min(length(sig_src),length(sig_rec))));
[~,index] = max(Cm);
index = index+1;

%% Trim Signals according to Xcorr results
if length(sig_rec) < length(sig_src)
%     sig_rec = sig_rec(index:end);
%     sig_src = sig_src(1:length(sig_rec));
%     sig_src = sig_src(index:end);
%     sig_rec = sig_rec(1:length(sig_src));
    sig_src = sig_src(index:min(length(sig_src),length(sig_rec)));
    sig_rec = sig_rec(1:min(length(sig_src),length(sig_rec)));
else
%     sig_src = sig_src(index:end);
%     sig_rec = sig_rec(1:length(sig_src));
    sig_rec = sig_rec(index:min(length(sig_rec),length(sig_src)));
    sig_src = sig_src(1:min(length(sig_rec),length(sig_src)));
end

%% Phase offset correction
phase_offset = angle(mean(sig_rec.*conj(sig_src)));
sig_rec = sig_rec*exp(-1i*phase_offset);

%% Fractional Delay Correction
for tt = 1:3
    [C1,~] = xcorr(sig_rec,sig_src,1);
    R = abs(C1);
    a = (R(3)+R(1))/2-R(2);
    b = (R(3)-R(1))/2;
    c = R(2);
    frac_delay = -b/(2*a);
    sig_rec = DelaySignalByFraction(sig_rec.',-frac_delay).';
end
%% sig_out Allocation
sig_out = sig_rec;
sig_out = sig_out/rms(sig_out)*rms(sig_src);
end
%% Auxilary Functions
function [DelayedSignal] = DelaySignalByFraction(Signal , Fraction)
n=length(Signal);
Delta = 2*pi/n;
if (mod(n,2)==1)
    w = [0  Delta:Delta:pi  fliplr(-Delta:-Delta:-pi)];
else
    w = [0:Delta:pi-Delta  fliplr(-Delta:-Delta:-pi)];
end
Phi = Fraction*w;
DelayedSignal = ifft(fft(Signal) .* exp(-1i*Phi));
end

