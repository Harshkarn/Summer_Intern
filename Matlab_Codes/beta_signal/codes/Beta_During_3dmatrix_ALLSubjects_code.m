clc
clear all;
folder_path='C:\Users\hkarn\OneDrive\Desktop\Data_during';
file_list=dir(fullfile(folder_path,'*.edf'));
During=zeros(31000,19,length(file_list));
for m=1:length(file_list)
    disp("Processing for Subject:")
    disp(m)
    file_name=file_list(m).name;
    full_file_path=fullfile(folder_path,file_name);
    [hdr272, record272] = edfread(full_file_path);
    %reads data from patient 27_2
beta_test=zeros(31000,19);
% Plotting all channels on separate subplots
for ch = 1:19
    
fz_272=record272(ch,:);% check channel fz
records_272=hdr272.records;
duration_272=hdr272.duration; %1sec

fs = 500;   %frequency
t=0:1/fs:records_272*duration_272-1/fs; %recording time

fz_272_fft=fft(fz_272); 
f_axis=0:fs/(length(fz_272)-1):fs/2; 

% figure(1)
% subplot(211)
% plot(t,fz_272);title('272 time domain (Channel fz)');grid on;
% xlabel('time [sec]');
% ylabel('amplitude [volt]');
% 
% subplot(212);plot(f_axis,fz_272_fft(1:length(fz_272_fft)/2));grid on;title('272 frequency domain (Channel fz)');
% xlabel('frequency[Hz]');
% ylabel('Amplitude[Volt]');

% section 3,4,5 --> Perform a wavelet decomposition 

[c272,l272] = wavedec(fz_272,6,'db8');

%detail coefficients
[cd1,cd2,cd3,cd4,cd5,cd6] = detcoef(c272,l272,[1 2 3 4 5 6]); 

% app coefficients

ca1 = appcoef(c272,l272,'db8',1);
ca2 = appcoef(c272,l272,'db8',2);
ca3 = appcoef(c272,l272,'db8',3);
ca4 = appcoef(c272,l272,'db8',4);
ca5 = appcoef(c272,l272,'db8',5);
ca6 = appcoef(c272,l272,'db8',6);
 
A1=wrcoef('a',c272,l272,'db8',1);
A2=wrcoef('a',c272,l272,'db8',2);
A3=wrcoef('a',c272,l272,'db8',3);
A4=wrcoef('a',c272,l272,'db8',4);
A5=wrcoef('a',c272,l272,'db8',5);
A6=wrcoef('a',c272,l272,'db8',6);


D1=wrcoef('d',c272,l272,'db8',1);
D2=wrcoef('d',c272,l272,'db8',2);
D3=wrcoef('d',c272,l272,'db8',3);
D4=wrcoef('d',c272,l272,'db8',4);
D5=wrcoef('d',c272,l272,'db8',5);
D6=wrcoef('d',c272,l272,'db8',6);


% section 7: time and freq domains for coefficients

%define time
t1=0:1/fs:(length(ca1)-1)/fs;
t2=0:1/fs:(length(ca2)-1)/fs;
t3=0:1/fs:(length(ca3)-1)/fs;
t4=0:1/fs:(length(ca4)-1)/fs;
t5=0:1/fs:(length(ca5)-1)/fs;
t6=0:1/fs:(length(ca6)-1)/fs;

ff_a1=0:fs/(2*(length(A1)-1)):fs/4;
Fca1=fft(A1);
ff_a2=0:fs/(4*(length(A2)-1)):fs/8;
Fca2=fft(A2);
ff_a3=0:fs/(8*(length(A3)-1)):fs/16;
Fca3=fft(A3);
ff_a4=0:fs/(16*(length(A4)-1)):fs/32;
Fca4=fft(A4);
ff_a5=0:fs/(32*(length(A5)-1)):fs/64;
Fca5=fft(A5);
ff_a6=0:fs/(64*(length(A6)-1)):fs/128;
Fca6=fft(A6);

ff_d1=fs/4:fs/(2*(length(D1)-1)):fs/2;
Fd1=fft(D1);
ff_d2=fs/8:fs/(4*(length(D2)-1)):fs/4;
Fd2=fft(D2);
ff_d3=fs/16:fs/(8*(length(D3)-1)):fs/8;
Fd3=fft(D3);
ff_d4=fs/32:fs/(16*(length(D4)-1)):fs/16;
Fd4=fft(D4);
ff_d5=fs/64:fs/(32*(length(D5)-1)):fs/32;
Fd5=fft(D5);
ff_d6=fs/128:fs/(64*(length(D6)-1)):fs/64;
Fd6=fft(D6);

% section 8:

[M,I]=max(Fd4(1:length(Fd4)/2));
fprintf('Beta:Maximum occurs at %3.2f Hz.\n',ff_d4(I));


beta1 = bandpass(D4,[13 30],500);
beta_test(:,ch)=beta1(1:31000)';

% Energy for Beta 1
j=1;
Ebeta1=0;
while (j<=length(beta1))
  Ebeta1=abs(beta1(j))^2+Ebeta1;
  j=j+1;
end



% section 10 

%reconstruct a noiseless signal

reconstSig=(D1+D2+D3+D4+D5+D6+A6); 
end
During(:,:,m)=beta_test;

% figure(5);
% plot(t,reconstSig);
% title('Noiseless Reconstructed Signal in Time');xlabel('time [sec]');ylabel('Amplitude [volt]')



% section 11 --> Denoising Original Signal

SIGDEN = cmddenoise(fz_272,'db8',6);
% plot(t,SIGDEN);title('SIGDEN in Time');xlabel('Time[sec]');ylabel('Amplitude[Volt]');
% grid on;


% Calculation of errors

err1=norm(SIGDEN-reconstSig); 
err2=norm(SIGDEN-fz_272); 
err3=norm(reconstSig-fz_272);
fprintf('the norm error between SIGDEN-reconstSig is %3.12f.\n',err1); 
fprintf('the norm error between SIGDEN-s is %3.12f.\n',err2); 
fprintf('the norm error between reconstSig-s is %3.12f.\n',err3);

end