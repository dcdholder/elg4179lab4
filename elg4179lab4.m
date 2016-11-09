function elg4179lab4()
	NUM_SNR_SAMPLES = 20;
	NUM_BIT_SAMPLES = 10000;

	snr = linspace(0,10,NUM_SNR_SAMPLES);
	
	%BPSK
	for n = 1:NUM_SNR_SAMPLES
		sampleData        = randomBitSource(NUM_SAMPLES);
		modulatedSignal   = modulator(sampleData);
		noisySignal       = noiseAddition(1/(10*log10(snr(n)),modulatedSignal);
		demodulatedOutput = demodulatorBPSK(noisySignal);
		errorRateBPSK(n)  = getErrorRate(sampleData,demodulatedOutput);
	end

	theoreticalErrorBPSK = 0.5*erfc(sqrt(1/(10*log10(snr)));
	
	figure(1);
	plot(snr,errorRateBPSK);
	title("SNR vs. Bit Error Rate for BPSK");
	xlabel("SNR (dB)");
	ylabel("Bit Error Rate");

	%OOK
	for n = 1:NUM_SNR_SAMPLES
		modulatedSignal   = randomBitSource(NUM_SAMPLES); %no need for modulation
		noisySignal       = noiseAddition(1/2/(10*log10(snr(n))),modulatedSignal);
		demodulatedOutput = demodulatorOOK(noisySignal);
		errorRateOOK(n)   = getErrorRate(modulatedSignal,demodulatedOutput); %"modulated" signal is a binary stream
	end

	theoreticalErrorOOK = 0.5*erfc(sqrt(1/2/(10*log10(snr)));

	figure(2);
	plot(snr,errorRateOOK);
	title("SNR vs. Bit Error Rate for OOK");
	xlabel("SNR (dB)");
	ylabel("Bit Error Rate");

	%QPSK
	for n = 1:NUM_SNR_SAMPLES
		sampleData              = randomSymbolSource(2,NUM_SAMPLES);
		modulatedSignal         = symbolModulatorQPSK(sampleData);
		noisySignal             = noiseAddition(1/(10*log10(snr(n)),modulatedSignal);
		demodulatedOutput       = demodulatorQPSK(noisySignal);
		errorRateQPSKSymbol(n)  = getSymbolErrorRateSymbol(2,sampleData,demodulatedOutput);
		errorRateQPSKBit(n)     = getBitErrorRateSymbol(2,sampleData,demodulatedOutput);
	end

	theoreticalErrorQPSKSymbol = erfc(sqrt(2*(10*log10(snr)));

	figure(3);
	plot(snr,errorRateQPSKSymbol);
	title("SNR vs. Symbol Error Rate for QPSK");
	xlabel("SNR (dB)");
	ylabel("Symbol Error Rate");
	
	theoreticalErrorQPSKBit = erfc(sqrt(4*(10*log10(snr)));
	
	figure(4);
	plot(snr,errorRateQPSKBit);
	title("SNR vs. Bit Error Rate for QPSK");
	xlabel("SNR (dB)");
	ylabel("Bit Error Rate");

	%16QAM
	for n = 1:NUM_SNR_SAMPLES
		sampleData              = randomSymbolSource(4,NUM_SAMPLES);
		modulatedSignal         = symbolModulator16QAM(sampleData);
		noisySignal             = noiseAddition(1/(10*log10(snr(n)),modulatedSignal);
		demodulatedOutput       = demodulator16QAM(noisySignal);
		errorRate16QAMSymbol(n) = getSymbolErrorRateSymbol(4,sampleData,demodulatedOutput);
		errorRate16QAMBit(n)    = getBitErrorRateSymbol(4,sampleData,demodulatedOutput);
	end

	theoreticalError16QAMSymbol = 3/2*erfc(sqrt(1/10*(10*log10(snr)));

	figure(5);
	plot(snr,errorRate16QAMSymbol);
	title("SNR vs. Symbol Error Rate for 16QAM");
	xlabel("SNR (dB)");
	ylabel("Symbol Error Rate");
	
	theoreticalError16QAMBit = 3/8*erfc(sqrt(8/5*(10*log10(snr)));
	
	figure(6);
	plot(snr,errorRate16QAMBit);
	title("SNR vs. Bit Error Rate for 16QAM");
	xlabel("SNR (dB)");
	ylabel("Bit Error Rate");
end

%random data sources

function bits = randomBitSource(numBits)
	bits = round(rand(1,numBits));
end

function symbols = randomSymbolSource(bitsPerSymbol,numSymbols)
	possibleSymbols = 2^bitsPerSymbol;
	symbols = floor(rand(possibleSymbols,numSymbols));
end

%modulators

function modulatedOutput = modulator(binaryData)
	modulatedOutput = binaryData.*2 -1;
end

function modulatedOutput = symbolModulatorQPSK(symbolData)
	for n=1:length(symbolData)
		switch symbolData(n)
			case 0:
				modulatedOutput(n) = -sqrt(1/2)-j*sqrt(1/2);
			case 1:
				modulatedOutput(n) = -sqrt(1/2)+j*sqrt(1/2);
			case 2:
				modulatedOutput(n) =  sqrt(1/2)-j*sqrt(1/2);
			case 3:
				modulatedOutput(n) =  sqrt(1/2)+j*sqrt(1/2);
		end
end

function modulatedOutput = symbolModulator16QAM(symbolData)
	for n=1:length(symbolData)
		switch symbolData(n)
			case 0:
				modulatedOutput(n) = sqrt(1/10)-j*sqrt(1/10);
			case 1:
				modulatedOutput(n) = 3*sqrt(1/10)-j*sqrt(1/10);
			case 2:
				modulatedOutput(n) = -sqrt(1/10)-j*sqrt(1/10);
			case 3:
				modulatedOutput(n) = -3*sqrt(1/10)-j*sqrt(1/10);
			case 4:
				modulatedOutput(n) = sqrt(1/10)-3*j*sqrt(1/10);
			case 5:
				modulatedOutput(n) = 3*sqrt(1/10)-3*j*sqrt(1/10);
			case 6:
				modulatedOutput(n) = -sqrt(1/10)-3*j*sqrt(1/10);
			case 7:
				modulatedOutput(n) = -3*sqrt(1/10)-3*j*sqrt(1/10);
			case 8:
				modulatedOutput(n) = sqrt(1/10)+j*sqrt(1/10);
			case 9:
				modulatedOutput(n) = 3*sqrt(1/10)+j*sqrt(1/10);
			case 10:
				modulatedOutput(n) = -sqrt(1/10)+j*sqrt(1/10);
			case 11:
				modulatedOutput(n) = -3*sqrt(1/10)+j*sqrt(1/10);
			case 12:
				modulatedOutput(n) = sqrt(1/10)+3*j*sqrt(1/10);
			case 13:
				modulatedOutput(n) = 3*sqrt(1/10)+3*j*sqrt(1/10);
			case 14:
				modulatedOutput(n) = -sqrt(1/10)+3*j*sqrt(1/10);
			case 15:
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
	demodulatedOutput = sign(real(modulatedSignalPlusNoise));
end

function demodulatedOutput = demodulatorOOK(modulatedSignalPlusNoise)
	demodulatedOutput = modulatedSignalPlusNoise-0.5;
	demodulatedOutput = sign(real(modulatedSignalPlusNoise));
end

function demodulatedSymbolOutput = demodulatorQPSK(modulatedSignalPlusNoise)
	for n=1:length(modulatedSignalPlusNoise)
		smallestDistance = 1000; %make this larger than the largest possible distance
		for compareDigit=0:3
			distance = abs(modulatedSignalPlusNoise(n)-symbolModulatorQPSK(compareDigit));
			if distance<smallestDistance
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
			numErrors = numErrors + bitsInError(numBits,originalInput(n),demodulatedOutput(n));
		end
	end
end

function numErrorBits = bitsInError(originalSymbol,compareSymbol)
	xorOut       = bitxor(originalSymbol,compareSymbol);
	maxSetBit    = ceil(log2(xorOut+0.1));
	numErrorBits = sum(bitget(xorOut,1:maxSetBit));
end
