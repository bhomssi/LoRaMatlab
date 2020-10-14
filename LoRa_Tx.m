function [signal_mod] = LoRa_Tx(message,Bandwidth,SF,Pt,Fs,df,varargin)
% LoRa_Tx emulates a Lora transmission
%
%   in:  message      payload message
%        Bandwidth    signal bandwidth of LoRa transmisson  
%        SF           spreading factor
%        Pt           transmit power in deicbels
%        Fs           sampling frequency
%        dF           frequency offset
%        varargin{1}  code rate
%        varargin{2}  symbols in preamble
%        varargin{3}  sync key
%
%  out:  signal       LoRa IQ waveform
%        packet       encoded message
%
% Dr Bassel Al Homssi  
% RMIT University 
% Credit to rpp0 on https://github.com/rpp0/gr-lora

if nargin == 6
    CR = 1 ;
    n_preamble = 8 ;
    SyncKey = 5 ;
elseif nargin == 7
    CR = varargin{1} ;
    n_preamble = 8 ;
    SyncKey = 5 ;
elseif nargin == 8
    CR = varargin{1} ;
    n_preamble = varargin{2} ;
    SyncKey = 5 ;
elseif nargin == 9
    CR = varargin{1} ;
    n_preamble = varargin{2} ;
    SyncKey = varargin{3} ;
end

packet = LoRa_Encode_Full(message,SF,CR) ; % encode message
signal = LoRa_Modulate_Full(packet,SF,Bandwidth,n_preamble,SyncKey,Fs) ; % LoRa modulate message
signal_mod = 10.^(Pt./20).*signal.*exp(-j.*2.*pi.*df/Fs.*(0:length(signal)-1))' ; % frquency shift and convert to power
end
function [packet] = LoRa_Encode_Full(message,SF,CR)
% LoRa_Encode_Full emulates a Lora transmission
%
%   in:  message      payload message
%        SF           spreading factor
%        CR           coding rate 
%
%  out:  packet       encoded lora packet 

CRC_pld = 1 ;  % cyclic rate code flag
imp = 0 ;
opt = 0 ;
%% String to Decimal
message_chr = convertStringsToChars(message) ;
message_dbl = uint8(message_chr) ;
%% Packet Length Calculations
N_pld = (SF == 7).*1 + (SF == 8).*2 + (SF == 9).*3 + (SF == 10).*4 + (SF == 11).*5 + (SF == 12).*6 ;
n_packet = 8 + max([ceil((8*(length(message_dbl) + 5) - 4.*SF + 28 + 16.*CRC_pld - 20.*imp)/(4.*(SF - 2.*opt))).*(CR + 4) 0]) ;
n_wht = SF .* floor((n_packet-8)/(4 + CR)) + N_pld - 1 ;
n_pld = ceil((n_wht + (SF == 7).*0 + (SF == 8).*1 + (SF == 9).*2 + (SF == 10).*3 + (SF == 11).*4 + (SF == 12).*5)/2) ;
n_pad = n_pld - 5 - length(message_dbl) - CRC_pld.*2 ;
%% Create payload message
CRC_dbl = CRC_pld.*[1 1] ; % CRC is not working atm
pad_dbl = zeros(1,n_pad + N_pld - 1) ; % padding
pld_dbl = [255 255 0 0 message_dbl 0 CRC_dbl pad_dbl] ; % LoRa payload
%% Swap Nibbles
pld_swp = LoRa_encode_swap(pld_dbl) ;
%% Payload Encoding
pld_enc = LoRa_encode_hamming(pld_swp,CR) ;
pld_enc = pld_enc(1 : n_wht) ;
%% Payload Whiten
pld_wht = LoRa_encode_white(pld_enc,CR,0) ;
%% Header Encoding
packet_hdr = [(length(message_dbl)+5) CRC_pld*16+(CR==1)*32+(CR==2).*64+(CR==3).*96+(CR==4)*128 224] ;
packet_hdr_enc_tmp = LoRa_encode_hamming(packet_hdr,4) ;
packet_hdr_enc = [packet_hdr_enc_tmp(1:5) pld_wht(1:N_pld-1)] ;
%% Packet Creation
packet_pld = pld_wht(N_pld : end) ;
packet_pld_shf = bitand(LoRa_encode_shuffle(packet_pld),2^(4+CR)-1)  ;
packet_hdr_shf = LoRa_encode_shuffle(packet_hdr_enc)  ;
% Interleaving
packet_pld_int = LoRa_encode_interleave(packet_pld_shf,SF,CR) ;
packet_hdr_int = LoRa_encode_interleave(packet_hdr_shf,SF-2,4) ;
% Graying
packet_pld_gray = LoRa_encode_gray(packet_pld_int) ;
packet_hdr_gray = LoRa_encode_gray(packet_hdr_int) ;
% Packet final
packet = [4*packet_hdr_gray packet_pld_gray] ;
end
function [symbols_swp] = LoRa_encode_swap(symbols)
% LoRa_encode_swap swaps nibbles
%
%   in:  symbols            symbol sequence
%
%  out:  symbols_swp        symbols with swapped nibbles

symbols_swp = zeros(1,length(symbols)) ;
for ctr = 1 : length(symbols)
    symbols_swp(ctr) = bitor(bitsll(bitand(symbols(ctr),hex2dec('0F')),4),bitsra(bitand(symbols(ctr),hex2dec('F0')),4)) ; % swap first half of 8-bit sequencne with other half 
end
end
function [encoded] = LoRa_encode_hamming(symbols,CR)
% LoRa_encode_hamming hamming encodes symbols to ensure a more accurate decoding
%
%   in:  symbols      symmbol sequence
%        CR           hamming coding rate 
%
%  out:  encoded      hamming encoded symbols

if CR > 2 && CR <= 4 % detection and correction
    n = floor(length(symbols).*(4 + 4)/4) ;
    
    H = [0,210,85,135,153,75,204,30,225,51,180,102,120,170,45,255] ; 
    
    encoded = zeros(1,n) ;
    Ctr = 1 ;
    for ctr = 1 : length(symbols)
        s0 = bitand(floor(bitsra(symbols(ctr),4)),hex2dec('0F')) ;
        s1 = bitand(floor(bitsra(symbols(ctr),0)),hex2dec('0F')) ;
        encoded(Ctr+0) = H(s0+1) ;
        encoded(Ctr+1) = H(s1+1) ;
        Ctr = Ctr + 2 ;
    end
elseif CR > 0 && CR <= 2 % detection
    Ctr = 1 ;
    for ctr = 1 : length(symbols)
        s0 = bitand(floor(bitsra(symbols(ctr),4)),hex2dec('FF')) ;
        s1 = bitand(floor(bitsra(symbols(ctr),0)),hex2dec('FF')) ;
        encoded(Ctr+0) = selectbits_encode(s0) ;
        encoded(Ctr+1) = selectbits_encode(s1) ;
        Ctr = Ctr + 2 ;
    end
end
end
function [symbols_white] = LoRa_encode_white(symbols,CR,DE)
% LoRa_encode_white symbols whitening by adding a known sequence to the payload
% bytes to reduce correlation redudancy
%
%   in:  symbols      symmbol sequence
%        CR           coding rate 
%        DE           data rate optimization flag
%
%  out:  symbols_white      whitened symbols

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
N = min([length(symbols) length(white_sequence)]) ; % cut-off to length of transmit symbols
symbols_white = bitxor(symbols(1:N),white_sequence(1:N)); % encode white
end
function [symbols_shuf] = LoRa_encode_shuffle(symbols)
% LoRa_encode_shuffle shuffles symbols by a to combine header and
% payload
%
%   in:  symbols            symbol vector
%
%  out:  symbols_shuf       shuffle symbols

for Ctr = 1 : length(symbols)
    symbols_binary = de2bi(symbols(Ctr),8) ;
    symbols_shuf_binary = [symbols_binary(2) symbols_binary(3) symbols_binary(4) ...
        symbols_binary(6) symbols_binary(5) symbols_binary(1) symbols_binary(7) ...
        symbols_binary(8)] ;
    symbols_shuf(Ctr) = bi2de(symbols_shuf_binary) ;
end
end
function [symbols_interleaved] = LoRa_encode_interleave(symbols,ppm,rdd)
% LoRa_encode_interleave imposes transposition and digit shift on the
% symbols and rotation
%
%   in:  symbols            symbol sequence
%        ppm                SF
%        rdd                CR
%
%  out:  symbols_interleaved       interleaved symbols

symbols_interleaved = [] ;
sym_idx_ext = 1 ;
for block_idx = 1 : floor(length(symbols)/(ppm))
    x = symbols((block_idx-1).*ppm+1:block_idx.*ppm) ;
    symbols_block_binary = de2bi(x,4+rdd) ;
    symbols_block_binary_rotated = transpose(symbols_block_binary) ; % transposed 
    symbols_block_rorated = bi2de(symbols_block_binary_rotated) ;
    mask = ppm ;
    % rotate
    for ctr = 1 : 4 + rdd
        sym_int(ctr) = rotl(symbols_block_rorated(ctr),mask,ppm) ;
        mask = mask - 1 ;
    end
    symbols_interleaved = [symbols_interleaved sym_int] ;
end
end
function [symbols] = LoRa_encode_gray(symbols)
% LoRa_encode_gray implements gray coding to reduce errors of adjacent bits.
%
%   in:  symbols            symbol sequence
%
%  out:  symbols            gray coded symbols

% XOR each symbol with a shifted mask
for ctr = 1 : length(symbols)
    symbols(ctr) = bitxor(symbols(ctr),floor(bitsra(symbols(ctr),16))) ;
    symbols(ctr) = bitxor(symbols(ctr),floor(bitsra(symbols(ctr),08))) ;
    symbols(ctr) = bitxor(symbols(ctr),floor(bitsra(symbols(ctr),04))) ;
    symbols(ctr) = bitxor(symbols(ctr),floor(bitsra(symbols(ctr),02))) ;
    symbols(ctr) = bitxor(symbols(ctr),floor(bitsra(symbols(ctr),01))) ;
end
end
function [symbol_rot] = selectbits_encode(symbol)
% selectbits_encode concat zeros (from 8-bit to 4-bit)
%
%   in:  symbols            symbol sequence
%
%  out:  symbols_rot        symbols for rotation

symbol_binary = de2bi(symbol,8) ;
symbol_binary_rot = [0 symbol_binary(1) symbol_binary(2) symbol_binary(3) 0 symbol_binary(4) 0 0] ;
symbol_rot = bi2de(symbol_binary_rot) ;
end
function [y] = rotl(bits,count,size)
% rotl modulo rotation
%
%   in:  bits            bit sequence
%        counts          
%        size
%
%  out:  y               symbols

len_mask = bitsll(1,size) - 1 ;
count = mod(count,size) ;
bits = bitand(bits,len_mask) ;
y = bitor(bitand(bitsll(bits,count),len_mask), floor(bitsra(bits,size - count))) ;
end
function [signal] = LoRa_Modulate_Full(packet,SF,Bandwidth,n_preamble,SyncKey,Fs)
% LoRa_Modulate_Full constructs a lora packet (preamble + sync header + payload)
%
%   in:  packet         payload 1xN symbol vector wher N=1-Inf 
%                       with values {0,1,2,...,2^(SF)-1}
%        SF             spreading factor   
%        Bandwidth      signal bandwidth of LoRa transmisson  
%        n_preamble     number of symbols in the preamble
%        SyncKey        synchronize key
%        Fs             sampling frequency
%
%  out:  signal          LoRa IQ packet

signal_prmb = loramod((SyncKey - 1).*ones(1,n_preamble),SF,Bandwidth,Fs,1) ; % preamble upchirps

signal_sync_u = loramod([0 0],SF,Bandwidth,Fs,1) ; % sync upchirp

signal_sync_d1 = loramod(0,SF,Bandwidth,Fs,-1) ; % header downchirp
signal_sync_d = [signal_sync_d1; signal_sync_d1; signal_sync_d1(1:length(signal_sync_d1)/4)] ; % concatenate header

signal_mesg = loramod(mod(packet + SyncKey,2^SF),SF,Bandwidth,Fs,1) ; % add sync key to payload messaage
signal = [signal_prmb; signal_sync_u; signal_sync_d; signal_mesg] ; % concatenate LoRa packet
end
function [y] = loramod(x,SF,BW,fs,varargin)
% loramod LoRa modulates a symbol vector specified by x
%
%   in:  x          1xN symbol vector wher N=1-Inf 
%                   with values {0,1,2,...,2^(SF)-1}
%        BW         signal bandwidth of LoRa transmisson  
%        SF         spreading factor   
%        Fs         sampling frequency
%        varargin{1} polarity of chirp
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
freq    = mod(gamma + Inv.*beta.*time,BW) - BW/2 ;

Theta   = cumtrapz(time,freq) ;
y       = reshape(exp(j.*2.*pi.*Theta),numel(Theta),1) ;
end
