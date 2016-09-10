
function main
clc;
close all;


%% FIR濾波算法

% 加載音頻
[y,fs,bits]=wavread('test');
% 獲取音頻長度
n=length(y);
% 頻譜變換
Y=fft(y);
% 繪製圖形
figure
subplot(2,1,1);plot(y);title('The original signal waveform');grid;
subplot(2,1,2);plot(abs(Y));title('Spectrum of the original signal');grid;
pause(2);
% 播放原始音頻
sound(y);
% 加入隨機噪聲
Noise=0.05*randn(n,1);
s=y+Noise;
S=fft(s);
% 繪製圖形
figure
subplot(2,1,1);plot(s);title('signal waveform(After adding noise)');grid;
subplot(2,1,2);plot(abs(S));title('signal spectrum(After adding noise)');grid;
pause(2);
% 播放加噪音頻
sound(s);

% FIR對音頻去噪
Fp1 = 800;
Fs1 = 1200;
Ft = fs;
wp=2*pi*Fp1/Ft;%通帶截止頻率 
ws = 2*pi*Fs1/Ft;%阻帶截止頻率
rp = 1;
rs = 50;
p = 1-10.^(-rp/20); %通帶阻帶波紋
s1 = 10.^(-rs/20);
fpts = [wp ws];
mag = [1 0];
dev = [p s1];
[n21,wn21,beta,ftype] = kaiserord(fpts,mag,dev);%kaiserord求階數截止頻率
b21 = fir1(n21,wn21,kaiser(n21+1,beta)); %由fir1設計濾波器
c = fftfilt(b21,s);
C=fft(c);
figure
subplot(2,1,1);plot(c);title('signal waveform(After FIR filter noise)');grid;
subplot(2,1,2);plot(abs(C));title('signal spectrum(After FIR filter noise)');grid;
pause(2);
% 播放FIR濾波後音頻
sound(c);
% 保存加噪音頻
wavwrite(s,fs,bits,'test+Noise.wav');
% 保存去噪音頻
wavwrite(c,fs,bits,'test+Noise+filter_FIR.wav');


%% IIR巴特沃斯濾波 
% 設置IIR濾波參數（截止頻率，阻帶頻率）  Butter可用在低通、高通、帶通和帶阻和模擬IIR濾波器
fc = 800;
fp = 1200;
Wp = 2*fc/fs;
Ws = 2*fp/fs;
Ap = 1; % 最小增益
As = 30;% 最大衰減
[N,wn]= buttord(Wp,Ws,Ap,As); % 設計濾波器 其中n代表濾波器階數,
%Wn代表波器的截止頻率,Wp和Ws分別是通帶和阻帶的截止頻率
%Rp和Rs分別是通帶和阻帶的波文關系
[b1,a1] = butter(N,wn); % 用巴特濾波器濾波
y_noise_filter = filter(b1,a1,s); % 用IIR濾波器濾波
C1 = fft(y_noise_filter); % 頻譜變換
figure
subplot(2,1,1);plot(y_noise_filter);title('signal waveform(After IIR filter noise)');grid;
subplot(2,1,2);plot(abs(C1));title('signal spectrum(After IIR filter noise)');grid;
pause(2);
% 播放IIR濾波後音頻
sound(y_noise_filter);
% 保存去噪音頻
wavwrite(y_noise_filter,fs,bits,'test+Noise+filter_IIR.wav');

%% 維納濾波
Rz=xcorr(s);
Gz=fft(Rz,n);
Rsz=xcorr(s,y);
Gsz=fft(Rsz,n);
Py = fftn(y);
H=Gsz./Gz; %維納濾波器的傳遞函數
S=H.*Py;

ss=real(ifft(S)); %原始信號的估計
ss=ss(1:n);
C2 = fft(ss);
figure
subplot(2,1,1);plot(ss);title('signal waveform(After weina filter noise)');grid;
subplot(2,1,2);plot(abs(C2));title('signal spectrum(After weina filter noise)');grid;
pause(2);
% 播放IIR濾波後音頻
sound(ss);
% 保存去噪音頻
wavwrite(y_noise_filter,fs,bits,'test+Noise+filter_weina.wav');


%% 小波算法濾波
%採用ddencmp函數獲得信號的默認閾值
[thr,sorth,keepapp]=ddencmp('den','wv',s);

%ddencmp結果顯示
thr
sorth
keepapp

% ddencmp為實現matlab中閾值獲取的函數之一，其的調用格式有以下三種：
% 1.[THR,SORH,KEEPAPP,CRIT]=ddencmp(IN1,IN2,X)
% 2.[THR,SORH,KEEPAPP,CRIT]=ddencmp(IN1,'wp',X)
% 3.[THR,SORH,KEEPAPP,CRIT]=ddencmp(IN1,'wv',X)
% 函數ddencmp用於獲取信號在消噪或壓縮過程中的默認閾值。
% 輸入參數X為一維或二維信號；
% IN1取值為'den'或'cmp'，'den'表示進行去噪，'cmp'表示進行壓縮；
% IN2取值為'wv'或'wp'，wv表示選擇小波，wp表示選擇小波包。
% 返回值THR是返回的閾值；SORH是軟閾值或硬閾值選擇參數；
% KEEPAPP表示保存低頻信號.

%%%%%採用wdencmp函數進行信號的閾值量化（即信號的消噪）%%
data1=wdencmp('gbl',s,'db5',5,thr,sorth,keepapp);

% 函數wdencmp的調用格式有以下三種：
% (1)[XC,CXC,LXC,PERF0,PERFL2]=wdencmp('gbl',X,'wname',N,THTR,SORH,KEEPAPP);
% (2)[XC,CXC,LXC,PERF0,PERFL2]=wdencmp('lvd',X,'wname',N,THTR,SORH);
% (3)[XC,CXC,LXC,PERF0,PERFL2]=wdencmp('lvd',C,L,'wname',N,THTR,SORH);
% 函數wdencmp用於一維或二維信號的消噪或壓縮。 wname是所用的小波函數，
% gbl(global的縮寫)表示每一層都採用同一個閾值進行處理，
% lvd表示每層採用不同的閾值進行處理，N表示小波分解的層數，
% THR為閾值向量，對於格式（2）和（3）每層都要求有一個閾值，因此閾值向量THR的長度為N，
% SORH表示選擇軟閾值或硬閾值（分別取值為's'和'h'），
% 參數KEEPAPP取值為1時，則低頻係數不進行閾值量化，反之，低頻係數要進行閾值量化。
% XC是要進行消噪或壓縮的信號，[CXC,LXC]是XC的小波分解結構，PERF0和PERFL2是恢復或壓縮L^2的範數百分比。

%繪製消噪後的的信號
C3 = fft(data1); % 頻譜變換
figure
subplot(2,1,1);plot(data1);title('signal waveform(After xiaobo filter noise)');grid;
subplot(2,1,2);plot(abs(C3));title('signal spectrum(After xiaobo filter noise)');grid;
pause(2);
% 播放小波濾波後音頻
sound(data1);
% 保存去噪音頻
wavwrite(data1,fs,bits,'test+Noise+filter_xiaobo.wav');




