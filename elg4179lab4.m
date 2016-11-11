function elg4179lab4()
	NUM_SNR_SAMPLES = 10;
	NUM_BIT_SAMPLES = 10000;

	snrDB = linspace(10/NUM_SNR_SAMPLES,10,NUM_SNR_SAMPLES);
	snr   = 10.^(snrDB./10);
	
	%BPSK
	for n = 1:NUM_SNR_SAMPLES
		sampleData        = randomBitSource(NUM_BIT_SAMPLES);
		modulatedSignal   = modulator(sampleData);
		noisySignal       = noiseAddition(1/(snr(n)),modulatedSignal);
		demodulatedOutput = demodulatorBPSK(noisySignal);
		errorRateBPSK(n)  = getErrorRate(sampleData,demodulatedOutput);
	end

	for n = 1:NUM_SNR_SAMPLES
		theoreticalErrorBPSK(n) = 0.5*erfc(sqrt(snr(n)));
	end
	
	figure(1);
	semilogy(snrDB,errorRateBPSK,'r',snrDB,theoreticalErrorBPSK,'b+');
	title('SNR vs. Bit Error Rate for BPSK');
	xlabel('SNR (dB)');
	ylabel('Bit Error Rate');
	legend('Measured','Theoretical');

	%OOK
	for n = 1:NUM_SNR_SAMPLES
		modulatedSignal   = randomBitSource(NUM_BIT_SAMPLES); %no need for modulation
		noisySignal       = noiseAddition(1/2/(snr(n)),modulatedSignal);
		demodulatedOutput = demodulatorOOK(noisySignal);
		errorRateOOK(n)   = getErrorRate(modulatedSignal,demodulatedOutput); %"modulated" signal is a binary stream
	end

	for n = 1:NUM_SNR_SAMPLES
		theoreticalErrorOOK(n) = 0.5*erfc(sqrt(snr(n)/2));
	end

	figure(2);
	semilogy(snrDB,errorRateOOK,'r',snrDB,theoreticalErrorOOK,'b+');
	title('SNR vs. Bit Error Rate for OOK');
	xlabel('SNR (dB)');
	ylabel('Bit Error Rate');
	legend('Measured','Theoretical');

	%QPSK
	for n = 1:NUM_SNR_SAMPLES
		sampleData              = randomSymbolSource(2,round(NUM_BIT_SAMPLES/2));
		modulatedSignal         = symbolModulatorQPSK(sampleData);
		noisySignal             = noiseAddition(1/2/snr(n),modulatedSignal);
		demodulatedOutput       = demodulatorQPSK(noisySignal);
		errorRateQPSKSymbol(n)  = getSymbolErrorRateSymbol(sampleData,demodulatedOutput);
		errorRateQPSKBit(n)     = getBitErrorRateSymbol(2,sampleData,demodulatedOutput);
	end

	for n = 1:NUM_SNR_SAMPLES
		theoreticalErrorQPSKSymbol(n) = erfc(sqrt(snr(n)));
	end
	
	figure(3);
	semilogy(snrDB,errorRateQPSKSymbol,'r',snrDB,theoreticalErrorQPSKSymbol,'b+');
	title('SNR vs. Symbol Error Rate for QPSK');
	xlabel('SNR (dB)');
	ylabel('Symbol Error Rate');
	legend('Measured','Theoretical');
	
	for n = 1:NUM_SNR_SAMPLES
		theoreticalErrorQPSKBit(n) = 0.5*erfc(sqrt(snr(n)));
	end
	
	figure(4);
	semilogy(snrDB,errorRateQPSKBit,'r',snrDB,theoreticalErrorQPSKBit,'b+');
	title('SNR vs. Bit Error Rate for QPSK');
	xlabel('SNR (dB)');
	ylabel('Bit Error Rate');
	legend('Measured','Theoretical');

	%16QAM
	for n = 1:NUM_SNR_SAMPLES
		sampleData              = randomSymbolSource(4,round(NUM_BIT_SAMPLES/4));
		modulatedSignal         = symbolModulator16QAM(sampleData);
		noisySignal             = noiseAddition(1/4/(snr(n)),modulatedSignal);
		demodulatedOutput       = demodulator16QAM(noisySignal);
		errorRate16QAMSymbol(n) = getSymbolErrorRateSymbol(sampleData,demodulatedOutput);
		errorRate16QAMBit(n)    = getBitErrorRateSymbol(4,sampleData,demodulatedOutput);
	end

	for n = 1:NUM_SNR_SAMPLES
		theoreticalError16QAMSymbol(n) = (3/2)*erfc(sqrt(2/5*snr(n)));
	end

	figure(5);
	semilogy(snrDB,errorRate16QAMSymbol,'r',snrDB,theoreticalError16QAMSymbol,'b+');
	title('SNR vs. Symbol Error Rate for 16QAM');
	xlabel('SNR (dB)');
	ylabel('Symbol Error Rate');
	legend('Measured','Theoretical');
	
	for n = 1:NUM_SNR_SAMPLES
		theoreticalError16QAMBit(n) = (3/8)*erfc(sqrt((2/5)*snr(n)));
	end
	
	figure(6);
	semilogy(snrDB,errorRate16QAMBit,'r',snrDB,theoreticalError16QAMBit,'b+');
	title('SNR vs. Bit Error Rate for 16QAM')
	xlabel('SNR (dB)');
	ylabel('Bit Error Rate');
	legend('Measured','Theoretical');
end

%random data sources

function bits = randomBitSource(numBits)
	bits = round(rand(1,numBits));
end

function symbols = randomSymbolSource(bitsPerSymbol,numSymbols)
	possibleSymbols = 2^bitsPerSymbol;
	symbols = floor(rand(1,numSymbols).*possibleSymbols);
end

%modulators

function modulatedOutput = modulator(binaryData)
	modulatedOutput = binaryData.*2 -1;
end

function modulatedOutput = symbolModulatorQPSK(symbolData)
	for n=1:length(symbolData)
		switch symbolData(n)
			case 0
				modulatedOutput(n) = -sqrt(1/2)-j*sqrt(1/2);
			case 1
				modulatedOutput(n) = -sqrt(1/2)+j*sqrt(1/2);
			case 2
				modulatedOutput(n) =  sqrt(1/2)-j*sqrt(1/2);
			case 3
				modulatedOutput(n) =  sqrt(1/2)+j*sqrt(1/2);
		end
	end
end

function modulatedOutput = symbolModulator16QAM(symbolData)
	for n=1:length(symbolData)
		switch symbolData(n)
			case 0
				modulatedOutput(n) = sqrt(1/10)-j*sqrt(1/10);
			case 1
				modulatedOutput(n) = 3*sqrt(1/10)-j*sqrt(1/10);
			case 2
				modulatedOutput(n) = -sqrt(1/10)-j*sqrt(1/10);
			case 3
				modulatedOutput(n) = -3*sqrt(1/10)-j*sqrt(1/10);
			case 4
				modulatedOutput(n) = sqrt(1/10)-3*j*sqrt(1/10);
			case 5
				modulatedOutput(n) = 3*sqrt(1/10)-3*j*sqrt(1/10);
			case 6
				modulatedOutput(n) = -sqrt(1/10)-3*j*sqrt(1/10);
			case 7
				modulatedOutput(n) = -3*sqrt(1/10)-3*j*sqrt(1/10);
			case 8
				modulatedOutput(n) = sqrt(1/10)+j*sqrt(1/10);
			case 9
				modulatedOutput(n) = 3*sqrt(1/10)+j*sqrt(1/10);
			case 10
				modulatedOutput(n) = -sqrt(1/10)+j*sqrt(1/10);
			case 11
				modulatedOutput(n) = -3*sqrt(1/10)+j*sqrt(1/10);
			case 12
				modulatedOutput(n) = sqrt(1/10)+3*j*sqrt(1/10);
			case 13
				modulatedOutput(n) = 3*sqrt(1/10)+3*j*sqrt(1/10);
			case 14
				modulatedOutput(n) = -sqrt(1/10)+3*j*sqrt(1/10);
			case 15
				modulatedOutput(n) = -3*sqrt(1/10)+3*j*sqrt(1/10);
		end
	end
end

%random noise source

function modulatedSignalPlusNoise = noiseAddition(N0,modulatedSignal)
	for n = 1:length(modulatedSignal)
		noise(n) = sqrt(N0/2) * (randn()+j*randn());
	end
	
	modulatedSignalPlusNoise = noise + modulatedSignal;
end

%demodulators

function demodulatedOutput = demodulatorBPSK(modulatedSignalPlusNoise)
	demodulatedOutput = real(modulatedSignalPlusNoise)>0;
end

function demodulatedOutput = demodulatorOOK(modulatedSignalPlusNoise)
	demodulatedOutput = real(modulatedSignalPlusNoise)>0.5;
end

%find the minimum possible distance from a legitimate QAM codepoint, match to that codepoint
function demodulatedSymbolOutput = demodulatorQPSK(modulatedSignalPlusNoise)
	for n=1:length(modulatedSignalPlusNoise)
		smallestDistance = 1000; %make this larger than the largest possible distance
		for compareDigit=0:3
			distance = abs(modulatedSignalPlusNoise(n)-symbolModulatorQPSK(compareDigit));
			if distance<smallestDistance
				smallestDistance = distance;
				demodulatedDigit = compareDigit;
			end
		end
		
		demodulatedSymbolOutput(n) = demodulatedDigit;
	end
end

function demodulatedSymbolOutput = demodulator16QAM(modulatedSignalPlusNoise)
	for n=1:length(modulatedSignalPlusNoise)
		smallestDistance = 1000; %make this larger than the largest possible distance
		for compareDigit=0:15
			distance = abs(modulatedSignalPlusNoise(n)-symbolModulator16QAM(compareDigit));
			if distance<smallestDistance
				smallestDistance = distance;
				demodulatedDigit = compareDigit;
			end
		end
		
		demodulatedSymbolOutput(n) = demodulatedDigit;
	end
end

%count errors

function errorRate = getErrorRate(originalInput,demodulatedOutput) %this only works for binary input
	sumInputOutput = originalInput + demodulatedOutput;

	numErrors = 0;
	for n=1:length(demodulatedOutput)
		if sumInputOutput(n)==1
			numErrors = numErrors+1;
		end
	end
	
	errorRate = numErrors / length(demodulatedOutput);
end

function errorRate = getSymbolErrorRateSymbol(originalInput,demodulatedOutput)
	numErrors = 0;
	for n=1:length(demodulatedOutput)
		if originalInput(n)!=demodulatedOutput(n)
			numErrors = numErrors+1;
		end
	end
	
	errorRate = numErrors / length(demodulatedOutput);
end

function errorRate = getBitErrorRateSymbol(numBits,originalInput,demodulatedOutput)
	numErrors = 0;
	for n=1:length(demodulatedOutput)
		if originalInput(n)!=demodulatedOutput(n)
			numErrors = numErrors + bitsInError(originalInput(n),demodulatedOutput(n));
		end
	end
	errorRate = numErrors/(length(originalInput)*numBits);
end

%xor mask sets all bits that are different between the two symbols
%we can then count the set bits
function numErrorBits = bitsInError(originalSymbol,compareSymbol)
	xorOut       = bitxor(originalSymbol,compareSymbol);
	if xorOut!=0
		maxSetBit    = ceil(log2(xorOut+0.1));
		numErrorBits = sum(bitget(xorOut,1:maxSetBit));
	else
		numErrorBits = 0;
	end
end
