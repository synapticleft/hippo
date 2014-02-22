function save_model(fname,m,p)

% Ahist=cat(3,Ahist,A); % Save a history of A?
fprintf(['\nWriting file: ' fname '...']);
if p.use_gpu
    m.A = double(m.A);
    m.D = double(m.D);
end
%eval(['save /media/Expansion Drive/redwood/cadieu-twolayer/state/' fname ' m p']);
save(['/media/work/hippocampus/state/' fname],'m','p');
fprintf(' Done.\n');
