function [newCurr,newVolt, newDcurr]=signalPrep(n,structInUse)
%     used to clean the noise in the current of the measurement by
%       removing frequency bins of insignifant values and reconstructing
%       the signal. also looks for and fixes human errors in bumps in
%       measurements

    current = structInUse.c(:,n);
    voltage = structInUse.v(:,n);
    fs = structInUse.fs;
    freq=structInUse.freq(n);
    meanDiff=mean(current(5*fs+1:8*fs))-mean(current(1:fs));
    
    if meanDiff>1e-4
        current=current(structInUse.errOffset+1:end);
        voltage=voltage(structInUse.errOffset+1:end);
    end
    
    
    Dcurr=getFourierDerivative(current,freq,fs);
    current=cleanNoise(current,freq,fs);
    
    newCurr=singleLoopAverage(current,freq,fs);
    newVolt=singleLoopAverage(voltage,freq,fs);
    newDcurr=singleLoopAverage(Dcurr,freq,fs);

end


function signal=singleLoopAverage(signal,freq,fs)
%     try
    signal=resample(signal,freq,1);
%     catch
%         freq
%     end
    tempSignal=zeros(fs,length(signal)/fs-2);
    for n=1+freq/10:length(signal)/fs-freq/10
        tempSignal(:,n-1)=signal((n-1)*fs+1:n*fs);
    end
    signal=sum(tempSignal,2)/size(tempSignal,2);

end

function newsignal = cleanNoise(signal,freq,fs)
%     get the clean frequency domain of signal
    fnewsignal=cleanInFourier(signal,freq,fs);
%     convert to real time domain
    newsignal=real(ifft(fnewsignal));
end

function derivative = getFourierDerivative(signal,freq,fs)
%   get the clean frequency domain of signal
    fnewsignal=cleanInFourier(signal,freq,fs);
%     simplify format for readability
    N=length(signal);
%     frequency vector
    w=(2*pi*[0:N/2-1, 0, -N/2+1:-1])';
%     take the derivative through jw
    dfnewsignal=1i*w.*fnewsignal;
%     convert to real time domain
    derivative=real(ifft(dfnewsignal));
end

function fnewsignal = cleanInFourier(signal,freq,fs)
%     get fft of signal
    fsignal=fft(signal);
    N=length(signal);
    nn=N/fs; % signal length compared to sampling freq.
%     dc=w(1), fundamental is w(nn*freq+1) and harmonics are w(nn*n*freq+1)
%   create storage for all harmonic data and dc data
%     find up to 10 harmonics, while within fs/2.56, the nyquist frequency.
    n=1;
    while n<=11 && (n-1)*freq<8192/2.56 
%         save indexes for harmonics, and for the dc, 1.
        indSave(n)=(n-1)*nn*freq+1;
        n=n+1;
    end

    L=n-1;
% right side of spectrum, left if fftshift is used, or accurate frequency
% vector with plotting. Same number of harmonics, hence length(indsave)
    for n=1:L-1
        indSave(n+L)=N-n*nn*freq+1;
    end
    
%     create a new signal with only zeros, no need to loop through
%     everything
    fnewsignal=zeros(size(fsignal));
%     loop through all indSave points
    for n=indSave
%         copy the data for given frequency bin to newsignal, only dc and
%         harmonics
        fnewsignal(n)=fsignal(n);
    end
end