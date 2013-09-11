function a = rotatePhase(a,theta)

if ~exist('theta','var')
    theta = -1*circ_mean(angle(a(:)))-2;
end

a = a.*exp(1i*theta);