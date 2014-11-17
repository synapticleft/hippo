% mexica() - Perform Independent Component Analysis (ICA) decomposition
%            of psychophysiological data using the ICA algorithm of 
%            Bell & Sejnowski (1995) with the natural gradient feature 
%            of Amari, Cichocki & Yang, or the extended-ICA algorithm 
%            of Lee, Girolami & Sejnowski (in press). Binary version of
%            runica() by Sigurd Enghoff (enghoff@salk.edu).
% Usage:
%      simply >> [weights,sphere] = mexica(data);
%       or    
%        else >> [weights,sphere,activations,bias,signs,lrates] ...
%                                 = mexica(data,'Key1',Value1',...);
% Input_Variable:
%
% data        = input data (chans,frames*epochs). 
%               Note: If data consists of multiple discontinuous epochs, 
%               each epoch should be separately baseline-zero'd using:
%                  >> data = rmbase(data,frames,basevector);
%
% NOTE (under DEC): Exiting from mexica() using Cnrtl-C fails to release memory.
%
% Optional_Keywords       Keyword_Values                  Default_Values
%
% 'pca'       = decompose a principal component subspace(default -> 0=off)
%               of the data. Value is the number of PCs to retain.
% 'sphering'  = flag sphering of data ('on'/'off')      (default -> 'on')
% 'weights'   = initial weight matrix                   (default -> eye())
%                    (Note: if 'sphering' 'off',        (default -> spher())
% 'lrate'     = initial ICA learning rate << 1          (default -> heuristic)
% 'block'     = ICA block size (integer << datalength)  (default -> heuristic)
% 'anneal'    = annealing constant (0,1] (defaults -> 0.90, or 0.98, extended)
%                         controls speed of convergence
% 'annealdeg' = degrees weight change for annealing     (default -> 70)
% 'stop'      = stop training when weight-change < this (default -> 1e-6)
% 'maxsteps'  = maximum number of ICA training steps    (default -> 512)
% 'bias'      = perform bias adjustment ('on'/'off')    (default -> 'on')
% 'momentum'  = momentum gain [0,1]                     (default -> 0)
% 'extended'  = [N] perform tanh() "extended-ICA" with kurtosis estimation 
%               every N training blocks. If N < 0, fix number of sub-Gaussian
%               components to -N [faster than N>0]      (default -> off)
% 'posact'    = make all component activations net-positive(default 'on'}
% 'verbose'   = give ascii messages ('on'/'off')        (default -> 'on')
%
% Output_Variables [RO = output in reverse order of projected mean variance 
%                        unless starting weight matrix passed ('weights' above)]
%
% weights     = ICA weight matrix (comps,chans)     [RO]
% sphere      = data sphering matrix (chans,chans) = spher(data)
%               Note: unmixing_matrix = weights*sphere {sphering off -> eye(chans)}
% activations = activation time courses of the output components 
%               (ncomps,frames*epochs)
% bias        = vector of final (ncomps) online bias [RO]    (default = zeros())
% signs       = extended-ICA signs for components    [RO]    (default = ones())
%                   [-1 = sub-Gaussian; 1 = super-Gaussian]
% lrates      = last learning rate used 

% Mex version:  Sigurd Enghoff, CNL / Salk Institute, La Jolla CA 7/98
%
% Toolbox Citation:
%
% Makeig, Scott et al. "ICA Toolbox for Psychophysiological Research (version 3.2)". 
% WWW Site, Computational Neurobiology Laboratory, The Salk Institute for Biological 
% Studies <www.cnl.salk.edu/~ica.html>, 1998. [World Wide Web Publication]. 
%
% Best Paper Reference:
%
% Scott Makeig, T-P. Jung, A. J. Bell, D. Ghahremani, T. J. Sejnowski,
% "Blind separation of auditory event-related brain responses into independent 
% components," Proceedings of the National Academy of Sciences USA, 
% 94:10979-10984 (1997).
%
% First Reference:
%
% Makeig, S., Bell, A.J., Jung, T-P and Sejnowski, T.J.,
% "Independent component analysis of electroencephalographic data," 
% In: D. Touretzky, M. Mozer and M. Hasselmo (Eds). Advances in Neural 
% Information Processing Systems 8:145-151, MIT Press, Cambridge, MA (1996).
%
% For more information:
% http://www.cnl.salk.edu/~scott/icafaq.html - FAQ on ICA/EEG
% http://www.cnl.salk.edu/~scott/icabib.html - mss. on ICA & biosignals
% http://www.cnl.salk.edu/~tony/ica.html - math. mss. on ICA, with kernal code
