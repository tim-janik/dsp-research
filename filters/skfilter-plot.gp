# lowpass 2p
HLP(s,k)=1/(s*s+2*(1-k/2)*s+1)

# bandpass 2p
HBP(s,k)=s/(s*s+2*(1-k/2)*s+1)

# highpass 2p
HHP(s,k)=s*s/(s*s+2*(1-k/2)*s+1)

# lowpass/highpass 1p
HHP1(s,k)=(s+1)*s/(s*s+2*(1-k/2)*s+1)
HLP1(s,k)=(s+1)/(s*s+2*(1-k/2)*s+1)

# plot [20:24000] db(abs(HLP(x/1000*{0,1.},0.95*2))), "<testfilter sweep 1000 4 0 0.95" w lines
