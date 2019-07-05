clear
clc
close all

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
d0 = 80;
v0 = -50;

%% FMCW Waveform Generation

%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
d_res = 1;
c = 3*10^8;
RMax = 200;

Bsweep = c/(2*d_res);  %Bandwidth
Tchirp = 5.5*2*RMax/c; %chirp time
alpha = Bsweep/Tchirp; %slope of chirps
fc= 77e9;              %carrier frequency of Radar 
                                                       
Nd=128;                %The number of chirps in one sequence
Nr=1024;               %The number of samples on each chirp.      

t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t));    %transmitted signal
Rx=zeros(1,length(t));    %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = d0 + v0*t(i);
    td(i) = 2*r_t(i)/c;
    
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*(fc*t(i)+alpha*t(i)^2/2));
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i)) + (alpha*(t(i)-td(i))^2)/2));
    
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);  %Beat Signal
    
end
%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.

%run the FFT on the beat signal along the range bins dimension (Nr) andnormalize.
sig_fft = fft(Mix,Nr)./Nr;

 % *%TODO* :
sig_fft = abs(sig_fft);       % Take the absolute value of FFT output

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
sig_fft = sig_fft(1:(Nr/2));

figure ('Name','Range from First FFT')  %plotting the range
plot(sig_fft); grid minor               % plot FFT output 
axis ([0 200 0 1]);
xlabel('measured range');

%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM

% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure('Name','Range Doppler Map'),surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tcr = 10;
Tcd = 4;

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gcr = 5;
Gcd = 2;

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 1.4;

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(Nr/2-2*(Tcd+Gcd),Nd-2*(Tcr+Gcr));
gridSize = (2*Tcr+2*Gcr+1)*(2*Tcd+2*Gcd+1);
trainingCellsNum = gridSize-(2*Gcr+1)*(2*Gcd+1);

% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.

CFAR_sig = zeros(size(RDM));

   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR
for j=1:Nd-2*(Tcr+Gcr)
    for i=1:Nr/2-2*(Tcd+Gcd)
        %I want to extract only the Training cells. To do that, I at first get
        %the sliding patch and, after converting it to power from decibel,
        %I set to zero whatever it is not in the position of
        %training cells; the zero submatrix will not contribute to the sum
        %over the whole patch and, by dividing for the number of training 
        %cells I will get the noise mean level
        
        trainingCellsPatch = db2pow(RDM(i:i+2*(Tcd+Gcd),j:j+2*(Gcr+Tcr)));
        trainingCellsPatch(Tcd+1:end-Tcd,Tcr+1:end-Tcr) = 0;
        
        noise_level(i,j) = pow2db(sum(sum(trainingCellsPatch))/trainingCellsNum);
        sigThresh = noise_level(i,j)*offset;
        if RDM(i+(Tcd+Gcd),j+(Tcd+Gcr))>sigThresh
            CFAR_sig(i+(Tcd+Gcd),j+(Tcd+Gcr)) = 1;
        else
            CFAR_sig(i+(Tcd+Gcd),j+(Tcd+Gcr)) = 0;
        end
           
    end
end


% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 


%not necessary



% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure('Name','CA-CFAR Filtered RDM'),surf(doppler_axis,range_axis,CFAR_sig);
colorbar;


 
 