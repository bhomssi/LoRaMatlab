function [message,symbols_Demod] = LoRa_Rx(signal,Bandwidth,SF,Coherece,Fs,df,varargin)
%LORA_RX Summary of this function goes here
%   Detailed explanation goes here
if nargin == 6
    SNR = Inf ;
elseif nargin == 7
    SNR = varargin{1} ;
end
if Fs == Bandwidth
    signal_demod = awgn(signal,SNR,'measured') ;
else
    signal_freq_demod = signal.*exp(j.*2.*pi.*df./Fs.*(0:length(signal)-1))' ;
    signal_filter = lowpass(signal_freq_demod,Bandwidth,Fs) ;
    signal_demod = awgn(resample(signal_filter,Bandwidth,Fs),SNR,'measured') ;
end
try
    [symbols_message,symbols_Demod] = LoRa_Demodulate_Full(signal_demod,SF,Bandwidth,Coherece);
    [message_full] = LoRa_Decode_Full(symbols_message,SF);
    message = message_full(8:4 + message_full(1) - 2) ;
catch
    message = NaN ;
end
end
function [symbols_message,symbols_Demod,n_preamble] = LoRa_Demodulate_Full(signal,SF,Bandwidth,Coherece)
%UNTITLED45 Summary of this function goes here
%   This code is to demodulate the signal
if SF > 12 || SF < 7
    return
end
n_symbol = 2^SF ;
%% Demodualte and Extract Preamble
n_signal = floor(length(signal)/n_symbol) ;
upChirps_demod = loramod(zeros(1,n_signal),SF,Bandwidth,Bandwidth) ;

% sniff_signal = signal(1:length(upChirps_demod)).*upChirps_demod ;
% fft_sync = zeros(1,n_signal) ;
% for Ctr = 1 : n_signal
%     fft_sync(Ctr) = FSKDetection(message,SF,Coherece) ;
% end
% [~,sync_ind] = sort(abs(fft_sync)) ;
% sync = sort(sync_ind(end-1:end)) ;
% sync = sync(end) + 1 ;
% n_preamble = sync - 5 ;
n_preamble = 8 ;
% disp(['Number of Preamble Symbols = ' num2str(n_preamble) ])

dnChirps_demod = loramod(zeros(1,n_preamble),SF,Bandwidth,Bandwidth,-1) ;
pream_signal = signal(1:length(dnChirps_demod)).*dnChirps_demod ;

symbols_pream = FSKDetection(pream_signal,SF,Coherece) ;
symbol_offset = mode(symbols_pream) + 1 ;
% disp(['Preamble Offset = ' num2str(symbol_offset) ])
%% Demodulate Message
message_start_ind = (n_preamble + 4.25)*n_symbol ;
n_message = length(signal)/n_symbol - message_start_ind/n_symbol ;
message_end_ind = n_message.*n_symbol + message_start_ind ;
% dnChirps_demod = loramod(zeros(1,n_message),SF,Bandwidth,Bandwidth,-1) ;
message = signal(message_start_ind+1:message_end_ind).*loramod(zeros(1,n_message),SF,Bandwidth,Bandwidth,-1) ;
symbols_Demod = FSKDetection(message,SF,Coherece) ;
symbols_message = mod(symbols_Demod - symbol_offset,2^SF) ;
% disp(['Message Symbols = ' num2str(symbols_message) ])
end
function [y] = loramod(x,SF,BW,fs,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if (nargin < 4)
    error(message('comm:pskmod:numarg1'));
end

if (nargin > 5)
    error(message('comm:pskmod:numarg2'));
end

% Check that x is a positive integer
if (~isreal(x) || any(any(ceil(x) ~= x)) || ~isnumeric(x))
    error(message('comm:pskmod:xreal1'));
end

M = 2^SF ;

% Check that M is a positive integer
if (~isreal(M) || ~isscalar(M) || M<=0 || (ceil(M)~=M) || ~isnumeric(M))
    error(message('comm:pskmod:Mreal'));
end

% Check that x is within range
if ((min(min(x)) < 0) || (max(max(x)) > (M-1)))
    error(message('comm:pskmod:xreal2'));
end

if nargin == 4
    Inv = 1 ;
elseif nargin == 5
    Inv = varargin{1} ;
end

Ts = 2^SF/BW ;
beta = BW/(2*Ts) ;
n_symbol = fs.*M/BW ;
t_symbol = (0:n_symbol-1).*1/fs ;

y = [] ;
for ctr = 1 : length(x)
    gamma = (x(ctr) - M/2)*BW/M ;
    lambda = 1 - x(ctr)/M ;
    t1 = t_symbol(1:end*lambda) ;
    t2 = t_symbol(end*lambda+1:end) ;
    y = [y; exp(-j.*2.*pi.*(t1'*gamma + beta*t1'.^2)*Inv); exp(-j.*2.*pi.*(t2'*(-BW + gamma) + beta*t2'.^2)*Inv)] ;
end
y = reshape(y,1,numel(y))' ;
end
function [message_full,CR_pld,pld_length,CRC_pld] = LoRa_Decode_Full(symbols_message,SF)
%LORA_DECODE Summary of this function goes here
%   Detailed explanation goes here
%% Decode Header
rdd_hdr = 4 ;
ppm_hdr = SF - 2 ;
symbols_hdr = mod(round(symbols_message(1:8)/4),2^ppm_hdr) ;
% Graying
symbols_hdr_gry = LoRa_decode_gray(symbols_hdr) ;
% Interleaving
symbols_hdr_int = LoRa_decode_interleave(symbols_hdr_gry,ppm_hdr,rdd_hdr) ;
% Shuffle
symbols_hdr_shf = LoRa_decode_shuffle(symbols_hdr_int,ppm_hdr) ;
% Hamming
symbols_hdr_fec = LoRa_decode_hamming(symbols_hdr_shf(1:5),rdd_hdr) ;

% disp(['Header Symbols = ' num2str(symbols_hdr_fec)])
%% Extract info from Header
CR_pld = floor(bitsra(symbols_hdr_fec(2),5)) ;
if CR_pld > 4 || CR_pld < 1
    return
end
CRC_pld = mod(floor(bitsra(symbols_hdr_fec(2),4)),2) ;
pld_length = symbols_hdr_fec(1) + CRC_pld*2 ;

% disp(['Payload Coding Rate = ' num2str(CR_pld)])
% disp(['Payload CRC = ' num2str(CRC_pld)])
% disp(['Payload Symbol Length = ' num2str(pld_length)])

% if SF > 10 || Data_Optimization == 1
%     flag = 2 ;
% else
%     flag = 0 ;
% end

% symbols_needed = floor((CR_pld + 4)*(8*pld_length)/(4*(SF - flag))) ;
% blocks_needed = ceil(symbols_needed/(CR_pld + 4)) ;
% pld_symbols = blocks_needed * (CR_pld + 4) ;
%% Decode Payload
rdd_pld = CR_pld ;

% if Data_Optimization == 1
%     ppm_pld = SF - 2 ;
%     symbols_pld = mod(round(symbols_message(9:end)/4),2^ppm_pld) ;
% else
ppm_pld = SF ;
symbols_pld = symbols_message(9:end) ;
% end

% Graying
symbols_pld_gry = LoRa_decode_gray(symbols_pld) ;
% Interleaving
symbols_pld_int = LoRa_decode_interleave(symbols_pld_gry,ppm_pld,rdd_pld) ;
% Shuffle
symbols_pld_shf = LoRa_decode_shuffle(symbols_pld_int,length(symbols_pld_int)) ;
% Add part of header
symbols_pld_hdr = [(SF>7).*symbols_hdr_shf(end - SF + 8:end) symbols_pld_shf] ;
% White
symbols_pld_wht = LoRa_decode_white(symbols_pld_hdr,rdd_pld,0) ;
% Hamming
symbols_pld_fec = LoRa_decode_hamming(symbols_pld_wht,rdd_pld) ;
% Swaping
symbols_pld_fin = LoRa_decode_swap(symbols_pld_fec) ;

% disp(['Payload Symbols = ' num2str(symbols_pld_swp(1:pld_length))])
%% Final Message
message_full = [symbols_hdr_fec symbols_pld_fin] ;
end
function [symbols_gray] = LoRa_decode_gray(symbols)
%LORA_DECODE_GRAY Summary of this function goes here
% XOR each symbol with a shifted mask
symbols_gray = bitxor(symbols,floor(bitsra(symbols,1))) ;
end
function [deocded] = LoRa_decode_hamming(symbols,CR)
%LORA_DECODE_HAMMING Summary of this function goes here
%   Detailed explanation goes here

if CR > 2 && CR <= 4
    n = ceil(length(symbols).*4/(4 + 4)) ;
    
    H = [0,0,0,0,0,0,3,3,0,0,5,5,14,14,7,7,0,0,9,9,2,2,7,7,4,4,7,7,7,7, ...
        7,7,0,0,9,9,14,14,11,11,14,14,13,13,14,14,14,14,9,9,9,9,10,10,9, ...
        9,12,12,9,9,14,14,7,7,0,0,5,5,2,2,11,11,5,5,5,5,6,6,5,5,2,2,1,1, ...
        2,2,2,2,12,12,5,5,2,2,7,7,8,8,11,11,11,11,11,11,12,12,5,5,14,14, ...
        11,11,12,12,9,9,2,2,11,11,12,12,12,12,12,12,15,15,0,0,3,3,3,3,3, ...
        3,4,4,13,13,6,6,3,3,4,4,1,1,10,10,3,3,4,4,4,4,4,4,7,7,8,8,13,13, ...
        10,10,3,3,13,13,13,13,14,14,13,13,10,10,9,9,10,10,10,10,4,4,13, ...
        13,10,10,15,15,8,8,1,1,6,6,3,3,6,6,5,5,6,6,6,6,1,1,1,1,2,2,1,1, ...
        4,4,1,1,6,6,15,15,8,8,8,8,8,8,11,11,8,8,13,13,6,6,15,15,8,8,1,1, ...
        10,10,15,15,12,12,15,15,15,15,15,15] ;
    
    deocded = zeros(1,n) ;
    for ctr = 0 : n - 1
        r0 = bitand(symbols(2*ctr+1),hex2dec("FF")) ;
        if 2*ctr+2 > length(symbols)
            symbols(2*ctr+2) = 0 ;
        end
        r1 = bitand(symbols(2*ctr+2),hex2dec("FF")) ;
        
        s0 = H(r0+1) ;
        s1 = H(r1+1) ;
        
        deocded(ctr+1) = bitor(bitsll(s0,4),s1) ;
    end
    
elseif CR > 0 && CR <= 2
    indices = [1 2 3 5] ;
    len = length(symbols) ;
    Ctr = 1 ;
    for ctr = 1 : 2 : len
        if ctr + 1 < len
            s1 = bitand(selectbits(symbols(ctr+1),indices),hex2dec("FF")) ;
        else
            s1 = 0 ;
        end
        s0 = bitand(selectbits(symbols(ctr),indices),hex2dec("FF")) ;
        deocded(Ctr) = bitor(bitsll(s0,4),s1) ;
        Ctr = Ctr + 1 ;
    end
end
end
function [symbols_interleaved] = LoRa_decode_interleave(symbols,ppm,rdd)
%UNTITLED42 Summary of this function goes here
%   Detailed explanation goes here
symbols_interleaved = [] ;
sym_idx_ext = 1 ;
for block_idx = 1 : floor(length(symbols)/(4+rdd))
    sym_int = zeros(1,ppm) ;
    for sym_idx = 1 : 4 + rdd
        sym_rot = rotl(symbols(sym_idx_ext),sym_idx-1,ppm) ;
        mask = bitsll(1,ppm-1) ;
        ctr = ppm ;
        while mask > 0
            sym_int(ctr) = sym_int(ctr) + bitsll(double(bitand(sym_rot,mask)>0),sym_idx-1) ;
            mask = floor(bitsra(mask,1)) ;
            ctr = ctr - 1 ;
        end
        sym_idx_ext = sym_idx_ext + 1 ;
    end
    symbols_interleaved = [symbols_interleaved sym_int] ;
end
end
function [symbols_shuf] = LoRa_decode_shuffle(symbols,N)
%LORA_DECODE_SHUFFLE Summary of this function goes here
%   Detailed explanation goes here
pattern = [5 0 1 2 4 3 6 7] ;
symbols_shuf = zeros(1,N) ;
for ctr = 1 : N
    for Ctr = 1 : length(pattern)
        symbols_shuf(ctr) = symbols_shuf(ctr) + bitsll(double(bitand(symbols(ctr),bitsll(1,pattern(Ctr)))>0),Ctr-1) ;
    end
end
end
function [symbols_swp] = LoRa_decode_swap(symbols)
%LORA_DECODE_SWAP Summary of this function goes here
%   Detailed explanation goes here
symbols_swp = zeros(1,length(symbols)) ;
for ctr = 1 : length(symbols)
    symbols_swp(ctr) = bitor(bitsll(bitand(symbols(ctr),hex2dec('0F')),4),bitsra(bitand(symbols(ctr),hex2dec('F0')),4)) ;
end
end
function [symbols_white] = LoRa_decode_white(symbols,CR,DE)
%LORA_DECODE_WHITE Summary of this function goes here
%   Detailed explanation goes here
if DE == 0
    if CR > 2 && CR <= 4
        white_sequence = [255,255,45,255,120,255,225,255,0,255,210,45,85, ...
            120,75,225,102,0,30,210,255,85,45,75,120,102,225,30,210,255, ...
            135,45,204,120,170,225,180,210,153,135,225,204,0,170,0,180,0, ...
            153,0,225,210,0,85,0,153,0,225,0,210,210,135,85,30,153,45,225, ...
            120,210,225,135,210,30,85,45,153,120,51,225,85,210,75,85,102, ...
            153,30,51,45,85,120,75,225,102,0,30,0,45,0,120,210,225,135,0, ...
            204,0,120,0,51,210,85,135,153,204,51,120,85,51,153,85,51,153, ...
            135,51,204,85,170,153,102,51,30,135,45,204,120,170,51,102,85, ...
            30,153,45,225,120,0,51,0,85,210,153,85,225,75,0,180,0,75,210, ...
            102,85,204,75,170,180,102,75,204,102,170,204,180,170,75,102, ...
            102,204,204,170,120,180,51,75,85,102,75,204,102,120,204,51, ...
            120,85,225,75,0,102,210,204,135,120,30,225,255,0,255,210,45, ...
            135,170,30,102,255,204,255,170,45,102,170,30,102,255,204,45, ...
            170,170,102,180,30,75,255,102,45,30,170,45,180,170,75,180,102, ...
            153,30,225,45,210,170,85,180,153,153,225,225,0,210,210,85,135, ...
            153,204,225,170,0,102,210,204,135,120,204,225,170,210,102,135, ...
            204,30,120,255,225,45,210,120,135,51,30,135,255,30,45,45,120, ...
            120,51,51,135,135,30,204,45,120,120,225,51,210,135,85,204,75, ...
            120,102,225,204,210,170,85,180,75,153,102,51,204,85,170,153, ...
            180,225,153,210,51,85,85,75,153,180,225,153,210,51,85,85,75, ...
            75,180,180,153,75,51,180,85,153,75,51,180,135,75,30,180,45, ...
            153,170,51,102,199,30,30,45,45,170,170,102,102,204,30,120, ...
            45,51,170,135,102,30,204,255,120,45,51,170,135,102,30,30,255, ...
            255,45,255,170,255,102,45,30,170,255,180,255,153,255,51,45,135, ...
            170,204,180,120,153,51,51,135,135,204,204,170,120,180,51,75,135, ...
            180,204,153,170,225,180,210,75,135,180,204,153,120,225,225,210,0, ...
            135,0,204,210,120,135,225,30,0,45,0,170,210,180,135,75,30,180, ...
            45,75,170,180,180,75,75,102,180,30,75,255,180,255,75,45,102,120, ...
            30,51,255,85,255,75,45,180,120,153,51,225,85,0,75,210,180,85,153, ...
            153,225,51,0,135,210,30,85,255,153,255,51,255,135,255,30,0,0,0,0, ...
            135,225,170,204] ;
    elseif CR > 0 && CR <= 2
        white_sequence = [255,255,45,255,120,255,48,46,0,46,18,60,20,40,10, ...
            48,54,0,30,18,46,20,60,10,40,54,48,30,18,46,6,60,12,40,58,48,36, ...
            18,24,6,48,12,0,58,0,36,0,24,0,48,18,0,20,0,24,0,48,0,18,18,6,20, ...
            30,24,60,48,40,18,48,6,18,30,20,60,24,40,34,48,20,18,10,20,54,24, ...
            30,34,60,20,40,10,48,54,0,30,0,60,0,40,18,48,6,0,12,0,40,0,34,18, ...
            20,6,24,12,34,40,20,34,24,20,34,24,6,34,12,20,58,24,54,34,30,6,60, ...
            12,40,58,34,54,20,30,24,60,48,40,0,34,0,20,18,24,20,48,10,0,36,0, ...
            10,18,54,20,12,10,58,36,54,10,12,54,58,12,36,58,10,54,54,12,12,58, ...
            40,36,34,10,20,54,10,12,54,40,12,34,40,20,48,10,0,54,18,12,6,40,30, ...
            48,46,0,46,18,60,6,58,30,54,46,12,46,58,60,54,58,30,54,46,12,60,58, ...
            58,54,36,30,10,46,54,60,30,58,60,36,58,10,36,54,24,30,48,60,18,58, ...
            20,36,24,24,48,48,0,18,18,20,6,24,12,48,58,0,54,18,12,6,40,12,48, ...
            58,18,54,6,12,30,40,46,48,60,18,40,6,34,30,6,46,30,60,60,40,40,34, ...
            34,6,6,30,12,60,40,40,48,34,18,6,20,12,10,40,54,48,12,18,58,20,36, ...
            10,24,54,34,12,20,58,24,36,48,24,18,34,20,20,10,24,36,48,24,18,34, ...
            20,20,10,10,36,36,24,10,34,36,20,24,10,34,36,6,10,30,36,60,24,58, ...
            34,54,6,30,30,60,60,58,58,54,54,12,30,40,60,34,58,6,54,30,12,46, ...
            40,60,34,58,6,54,30,30,46,46,60,46,58,46,54,60,30,58,46,36,46,24, ...
            46,34,60,6,58,12,36,40,24,34,34,6,6,12,12,58,40,36,34,10,6,36,12, ...
            24,58,48,36,18,10,6,36,12,24,40,48,48,18,0,6,0,12,18,40,6,48,30,0, ...
            60,0,58,18,36,6,10,30,36,60,10,58,36,36,10,10,54,36,30,10,46,36,46, ...
            10,60,54,40,30,34,46,20,46,10,60,36,40,24,34,48,20,0,10,18,36,20,24, ...
            24,48,34,0,6,18,30,20,46,24,46,34,46,6,46,30,0,0,0,0,36,6] ;
    end
end
N = min([length(symbols) length(white_sequence)]) ;
symbols_white = bitxor(symbols(1:N),white_sequence(1:N)) ;
end
function [y] = rotl(bits,count,size)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
len_mask = bitsll(1,size) - 1 ;
count = mod(count,size) ;
bits = bitand(bits,len_mask) ;
y = bitor(bitand(bitsll(bits,count),len_mask), floor(bitsra(bits,size - count))) ;
end
function [r] = selectbits(data,indices)
%SELECTBITS Summary of this function goes here
%   Detailed explanation goes here
r = 0 ;
for ctr = 0 : length(indices) - 1
    if bitand(data,bitsll(1,indices(ctr+1))) > 0
        r = r + bitsll(1,ctr) ;
    else
        r = r + 0 ;
    end
end
end
function [symbols] = FSKDetection(signal,SF,detection)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if detection == 1
    t = 0:1/(2^SF):0.999 ;
    for Ctr = 1 : 2^SF
        rtemp = conv(signal,exp(-j.*2.*pi.*(2^SF - Ctr + 1).*t)) ;
        r(Ctr,:) = real(rtemp(2^SF+1:2^SF:end)) ;
    end
    [~,idx] = max(r) ;
    symbols = idx - 1 ;
elseif detection == 2
    [~,idx] = max(fft(reshape(signal,2^SF,length(signal)/(2^SF)))) ;
    symbols = idx - 1 ;
end
end