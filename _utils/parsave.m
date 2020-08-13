function parsave(sname,S,varargin)
% parsave(sname,S,varargin)
% 
% wrapper fro save function, to run with parfor loops. "S" is a struct
% whose fields are the names of the variables to be saved

save(sname,'-struct','S',varargin{:})