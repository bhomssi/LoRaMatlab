function [message] = LoRa_Rx(signal,Bandwidth,SF,Coherece,Fs,df,varargin)
% LoRa_Rx emulates a Lora receiver
%
%   in:  signal       payload message
%        Bandwidth    signal bandwidth of LoRa transmisson  
%        SF           spreading factor
%        Coherence    (1) coherent or (2) non-coherent FSK Detection
%        Fs           sampling frequency
%        dF           carrier frequency shift
%        varagin{1}   SNR
%        varagin{2}   Preamble Symbol number
%
%  out:  message       LoRa payload message chahracters
%        symbols_Demod LoRa payload symbols vector  
%
% Dr Bassel Al Homssi  
% RMIT University 
% Credit to rpp0 on https://github.com/rpp0/gr-lora

if nargin == 6
    SNR                 = Inf ;
    n_preamble          = 8 ;
elseif nargin == 7
    SNR                 = varargin{1} ;
    n_preamble          = 8 ;
elseif nargin == 8
    SNR                 = varargin{1} ;
    n_preamble          = varargin{2} ;
end

if Fs == Bandwidth
    signal_demod        = awgn(signal,SNR,'measured') ;
else
    signal_freq_demod   = signal.*exp(j.*2.*pi.*df./Fs.*(0:length(signal)-1))' ;
    signal_filter       = lowpass(signal_freq_demod,Bandwidth,Fs) ;
    signal_demod        = awgn(resample(signal_filter,Bandwidth,Fs),SNR,'measured') ;
end

try
    symbols_message     = LoRa_Demodulate_Full(signal_demod,SF,Bandwidth,Coherece,n_preamble);
    [message_full]      = LoRa_Decode_Full(symbols_message,SF);
    message             = message_full(8:4 + message_full(1) - 2) ;
catch
    message             = NaN ;
end
end
function [SymbolsMessage,SymbolsDemod,NPreamb] = LoRa_Demodulate_Full(signal,SF,Bandwidth,Coherece,n_preamble)
% LoRa_Demodulate_Full demodulates full LoRa packet
%
%   in:  signal         IQ LoRa signal containing 
%                       (preamble + sync header + payload)
%        SF             spreading factor   
%        Bandwidth      signal bandwidth of LoRa transmisson  
%        Coherence      type of demodulation (coherent or non-coherent)
%
%  out:  symbols_message          message symbols vector (encoded)
%        symbols_Demod            LoRa symbols vector
%        n_preamble               Number of symbols in preamble
%% Return if SF is not in the range
if SF > 12 || SF < 7
    return
end
M = 2^SF ;
%% Demodualte and Extract Preamble
Nsymbols        = floor(length(signal)/M) ;
UChirpsDemod    = loramod(zeros(1,Nsymbols),SF,Bandwidth,Bandwidth) ;

SniffSignal     = signal(1:length(UChirpsDemod)).*UChirpsDemod ;

if Coherece == 2
    fftSync     = fft(reshape(SniffSignal,M,length(SniffSignal)/M)) ;
    [~,SyncInd] = sort(max(fftSync)) ;
    sync        = sort(SyncInd(end-1:end)) ;
    sync        = sync(end) + 1 ;
    NPreamb     = sync - 5 ;
else
    NPreamb     = n_preamble ;
end
dChirpsDemod    = loramod(zeros(1,NPreamb),SF,Bandwidth,Bandwidth,-1) ;
pream_signal    = signal(1:length(dChirpsDemod)).*dChirpsDemod ;

symbols_pream   = FSKDetection(pream_signal,SF,Coherece) ;
symbol_offset   = mode(symbols_pream) + 1 ;
%% Demodulate Message
MessageStartInd = (NPreamb + 4.25)*M ;
Nmessage        = floor(length(signal)/M - MessageStartInd/M) ;
MessageEndInd   = Nmessage.*M + MessageStartInd ;

MessageSignal   = signal(MessageStartInd+1:MessageEndInd).*loramod(zeros(1,Nmessage),SF,Bandwidth,Bandwidth,-1) ;
SymbolsDemod    = FSKDetection(MessageSignal,SF,Coherece) ;
SymbolsMessage  = mod(SymbolsDemod - symbol_offset,2^SF) ;
end
function [y] = loramod(x,SF,BW,fs,varargin)
% loramod LoRa modulates a symbol vector specified by x
%
%   in:  x          1xN symbol vector wher N=1-Inf 
%                   with values {0,1,2,...,2^(SF)-1}
%        BW         signal bandwidth of LoRa transmisson  
%        SF         spreading factor   
%        Fs         sampling frequency
%        varargin{1} set polarity of chirp
%
%  out:  y          LoRa IQ waveform
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

M       = 2^SF ;

% Check that M is a positive integer
if (~isreal(M) || ~isscalar(M) || M<=0 || (ceil(M)~=M) || ~isnumeric(M))
    error(message('comm:pskmod:Mreal'));
end

% Check that x is within range
if ((min(min(x)) < 0) || (max(max(x)) > (M-1)))
    error(message('comm:pskmod:xreal2'));
end

% Polarity of Chirp
if nargin == 4
    Inv = 1 ;
elseif nargin == 5
    Inv = varargin{1} ;
end
% Symbol Constants
Ts      = 2^SF/BW ;
Ns      = fs.*M/BW ;

gamma   = x/Ts ;
beta    = BW/Ts ;

time    = (0:Ns-1)'.*1/fs ;
freq    = mod(gamma + beta.*time,BW) - BW/2 ;

Theta   = cumtrapz(time,Inv.*freq) ;
y       = reshape(exp(j.*2.*pi.*Theta),numel(Theta),1) ;
end
function [message_full,CR_pld,pld_length,CRC_pld] = LoRa_Decode_Full(symbols_message,SF)
% LoRa_Decode_Full decodes full payload packet
%
%   in:  symbols_message         LoRa payload symbol vector
%        SF                      spreading factor  
%
%  out:  message_full          message symbols vector (decoded)
%        CR_pld                code rate of payload
%        pld_length            length of payload
%        CRC_pld               payload cyclic rate code flag
%% Decode Header
rdd_hdr         = 4 ;
ppm_hdr         = SF - 2 ;
symbols_hdr     = mod(round(symbols_message(1:8)/4),2^ppm_hdr) ;
% Graying
symbols_hdr_gry = LoRa_decode_gray(symbols_hdr) ;
% Interleaving
symbols_hdr_int = LoRa_decode_interleave(symbols_hdr_gry,ppm_hdr,rdd_hdr) ;
% Shuffle
symbols_hdr_shf = LoRa_decode_shuffle(symbols_hdr_int,ppm_hdr) ;
% Hamming
symbols_hdr_fec = LoRa_decode_hamming(symbols_hdr_shf(1:5),rdd_hdr) ;
%% Extract info from Header
CR_pld          = floor(bitsra(symbols_hdr_fec(2),5)) ;
if CR_pld > 4 || CR_pld < 1
    return
end
CRC_pld         = mod(floor(bitsra(symbols_hdr_fec(2),4)),2) ;
pld_length      = symbols_hdr_fec(1) + CRC_pld*2 ;
%% Decode Payload
rdd_pld         = CR_pld ;
ppm_pld         = SF ;
symbols_pld     = symbols_message(9:end) ;
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
%% Final Message
message_full    = [symbols_hdr_fec symbols_pld_fin] ;
end
function [symbols_gray] = LoRa_decode_gray(symbols)
% LoRa_decode_gray degray LoRa payload
%
%   in:  symbols       symbols with graying
%
%  out:  symbols_gray  degrayed symbols
symbols_gray = bitxor(symbols,floor(bitsra(symbols,1))) ;
end
function [deocded] = LoRa_decode_hamming(symbols,CR)
% LoRa_decode_hamming LoRa payload hamming decode (4,4 + CR)
%
%   in:  symbols       symbols with hamming
%        CR            Code Rate
%
%  out:  deocded      Fully decoded payload symbols

if CR > 2 && CR <= 4 % detection and correction
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
elseif CR > 0 && CR <= 2 % detection
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
% LoRa_decode_interleave deinterleaves payload packet
%
%   in:  symbols       interleaved symbols
%        ppm
%        rdd
%
%  out:  symbols_interleaved  deinterleaved symbols

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
% LoRa_decode_shuffle unshuffles payload packet
%
%   in:  symbols       symbol vector
%        N             
%
%  out:  symbols_shuf  unshuffled symbols

pattern = [5 0 1 2 4 3 6 7] ;
symbols_shuf = zeros(1,N) ;
for ctr = 1 : N
    for Ctr = 1 : length(pattern)
        symbols_shuf(ctr) = symbols_shuf(ctr) + bitsll(double(bitand(symbols(ctr),bitsll(1,pattern(Ctr)))>0),Ctr-1) ;
    end
end
end
function [symbols_swp] = LoRa_decode_swap(symbols)
% LoRa_decode_shuffle swap payload packet
%
%   in:  symbols       symbol vector           
%
%  out:  symbols_swp   unswapped symbols

symbols_swp = zeros(1,length(symbols)) ;
for ctr = 1 : length(symbols)
    symbols_swp(ctr) = bitor(bitsll(bitand(symbols(ctr),hex2dec('0F')),4),bitsra(bitand(symbols(ctr),hex2dec('F0')),4)) ; % swap first half of 8-bit sequencne with other half 
end
end
function [symbols_white] = LoRa_decode_white(symbols,CR,DE)
% LoRa_decode_white dewhitening of payload packet
%
%   in:  symbols       whitened symbols
%        CR            code rate
%        DE            data rate optimization flag
%
%  out:  symbols_white  dewhitened symbols

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
N = min([length(symbols) length(white_sequence)]) ; % LoRa symbol length
symbols_white = bitxor(symbols(1:N),white_sequence(1:N)) ;
end
function [y] = rotl(bits,count,size)
% rotl 
%
%   in:  bits            bit sequence
%        counts          
%        size
%
%  out:  y               rotated symbols

len_mask = bitsll(1,size) - 1 ;
count = mod(count,size) ;
bits = bitand(bits,len_mask) ;
y = bitor(bitand(bitsll(bits,count),len_mask), floor(bitsra(bits,size - count))) ;
end
function [r] = selectbits(data,indices)
% selectbits concat zeros (from 4-bit to 8-bit)
%
%   in:  data            symbol sequence
%        indices         vector = [1 2 3 4 5]
%
%  out:  r       `        symbols 

r = 0 ;
for ctr = 0 : length(indices) - 1
    if bitand(data,bitsll(1,indices(ctr+1))) > 0
        r = r + bitsll(1,ctr) ; % shift to left
    else
        r = r + 0 ;
    end
end
end
function [symbols] = FSKDetection(signal,SF,detection)
% LoRa_Tx demodulates a Lora de-chirped signal using
% the coherence specified by the detection variable
%
%   in:  message      payload message
%        SF           spreading factor
%        detection    1= coherent detection, 2= non-coherent detection   
%
%  out:  symbols      FSK demodulated symbol vector 

if detection == 1 % coherent detection
    t = 0:1/(2^SF):0.999 ; % time vector
    for Ctr = 1 : 2^SF
        rtemp = conv(signal,exp(-j.*2.*pi.*(2^SF - Ctr + 1).*t)) ; % convolution w/ideal fsk signal
        r(Ctr,:) = real(rtemp(2^SF+1:2^SF:end)) ; % save resultant array
    end
    [~,idx] = max(r) ; % take max
    symbols = idx - 1 ; % store symbol vector
elseif detection == 2 % non-coherent detection
    [~,idx] = max(fft(reshape(signal,2^SF,length(signal)/(2^SF)))) ; % take max of fft window
    symbols = idx - 1 ; % store symbol array
end
end