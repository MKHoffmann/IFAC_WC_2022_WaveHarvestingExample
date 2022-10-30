function [Value varargout] = HamiltonianPathForFBMNew(FBM,p)
arguments
    FBM     (:,:) {mustBeNumeric,mustBeCorrectDimension(FBM)}
    p       (1,1) {mustBeNumeric}
end
%% Returns the Value of the Simplex at the Position of the Path
% |O---->---->--|
% |xxxxxxxxxxxxx|
% ||----<----<--|
% |xxxxxxxxxxxxx|
% ||---->---->--|
% |xxxxxxxxxxxxx|
% ||----<----<--|
Vertex = size(FBM,1);   %Rectangular Simplex 
Pathlength = Vertex*Vertex-((((Vertex+1)/4)*2)-1)*(Vertex-2);
p = mod(p,Pathlength);

    
if p ==0
    p = Pathlength;
end
Row = ceil(p/(size(FBM,1)));
Row = (Row*2)-1;
if mod(p,size(FBM,1)) ==0
    Row = Row+1;
end


Column = mod(p,(size(FBM,1)));
if Column == 0
    [a b Column] = HamiltonianPathForFBMNew(FBM,p-1);
else 
    Column = Column+1;


end
if mod(Row+1,4)==0 
  
    Column = Vertex-Column+2;
    
end
if p > (Pathlength -Vertex);
Column = 1;
Row = Vertex-(p-(Pathlength -Vertex))+1;
end
if p ==1
    Row = 1;
    Column =2;
end
% For Debugging uncomment
% Simplex(Row,Column) = 0
Value = FBM(Row,Column);
varargout =cell(2);
varargout{1}=Row;
varargout{2} =Column;
end
% Custom validation function
function mustBeCorrectDimension(FBM)
    % Test for equal size
    if ~(mod(size(FBM,1)+1,4)==0)
        eid = 'Dimension:Wrong';
        msg = 'FBM Dimension must be: (N*4)-1';
        throwAsCaller(MException(eid,msg))
    end
end