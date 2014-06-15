ISEE-3/ICE Telemetry demodulator and decoder
Copyright Phil Karn, KA9Q, 15 June 2014
version 0.5
May be used under the terms of the GNU Public License 2.0 (GPL)


This is a snapshot of my ISEE-3/ICE telemetry decoder. It is now broken into four modules that you run
in a UNIX pipeline like this:

pmdemod input_file | symdemod | vdecode | framer

The input_file is expected to be a series of 16-bit signed integer samples in little-endian format. The I
channel is expected first. If Q is first, use the -f (flip) flag to invert the spectrum.

Each program writes its status to stderr so you can see what it's doing while its output is redirected.

This version has been tested at 32 sps/16 bps, 1024 sps/512 bps and 4096 sps/2048 bps. Note: it's important
to give a precise clock to symdemod with the -c option or symbol acquisition may not be reliable.
The measured clock speeds (at -3 km/s velocity) are:

1024 sps: 1024.475 Hz
4096 sps: 4097.9 Hz
32 sps: 32.0148 Hz

pmdemod - reads raw receiver file, tracks carrier, writes baseband PM samples (16-bit signed ints) on stdout

-S starting carrier frequency estimate, Hz
-W +/- Frequency search range while locked, Hz (searches entire passband when unlocked)
-D Doppler chirp rate, Hz/s
-t C/N0 lock threshold, dB; default 21 dB (good for 512 bps; should be lowered for lower data rates)
-q Quiet mode
-b carrier tracking FFT bin size, Hz; default 4 Hz (changed from 1 Hz to better track spin Doppler)
-f flip I & Q samples (invert spectrum)
-r Sample rate, Hz; default: 250000

symdemod - reads output of pmdemod on stdin, writes demodulated, 8-bit (offset-128) soft decisions on stdout
-c Symbol rate, Hz; default 1024.475 (measured rate at 1024 Hz nominal; divide or multiply by 2^n)
-C Clocks per symbol (default 1 for symbol rates > 1000; 1024/symrate for lower rates)
-r Sample rate, Hz; default 250000
-q quiet mode

vdecode - Reads output of symdemod on stdin, writes Viterbi-decoded bits on stdout ('0'/'1')
-p Start with opposite decoder symbol phase (also flips automatically based on encoded sync observations)
-i Status update interval, bits; default 1024
-d Viterbi decode delay, bits; default 200
-q quiet mode

framer - reads output of vdecode on stdin, detects frame sync, writes decoded telemetry frames in hex on stdout
-r bitrate, Hz; default 512 (only for time estimation)

The old programs icedemod, icesync and bitsync are included for reference but are not built by default.
