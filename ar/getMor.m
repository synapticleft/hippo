function w = getMor(freq,sd,rate,ns)
st = 1./(2*pi*sd);%df;%
w_sz = ns*st*rate; % half time window size
t = (-w_sz:(w_sz+1))/rate;
w = exp(-t.^2/(2*st^2)).*exp(2j*pi*freq*t)/sqrt(sqrt(pi)*st*rate);