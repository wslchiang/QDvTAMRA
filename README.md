# QDvTAMRA
Scripts and test files associated with QD fluorophore comparison paper

Summary of included files:

File Name						Format	Software to run		Description

Analysis_multi			.m			MATLAB						Analyze the QuickPALM outputs for simulations of multiple overlapped emitters
Analysis_single			.m			MATLAB						Analyzes the QuickPALM outputs for simulations of single emitters
blink_cp2		      	.m			MATLAB						Uses results from getData to perform changepoint trajectory and photoswitching analysis of traces
computeSNR		    	.m			MATLAB						Extracts the PL and BG values selected traces and computes the SNR for them
csv2png			      	.txt		Python						Plots localizations from QuickPALM that met a precision cutoff criterion. Points are plotted at center of mass with associated standard error
Extract_photophys		.mlx		MATLAB						Live script from MATLAB that calls upon other files to perform analysis of photophysical properties
getBkgSD		      	.m			MATLAB						Extracts the mean and standard deviation intensity values from BG traces
getData			      	.m			MATLAB						Extracts all necessary data points from PL and BG files and does initial processing of them to determine if trace has more than one emitter and computes the best intensity cutoff for "on" to "off" photoswitching
getPLSD			      	.m			MATLAB						Same as getBkgSD but for traces identified as PL
getResults					.m			MATLAB						Called upon by MakeFitHistograms to extract the desired values from the results of the analysis of getData and blink_cp2
MakeFitHistograms		.m			MATLAB						Plots the results for SNR and blink rate analysis
MultieEmitter				.m			MATALB						Generates series of images simulating multiple emitters with overlapped PSFs
SingleEmitter				.m			MATALB						Generates series of images simulating a single emitter with a noisy, blurred PSF
