function r = iscolinear(p1, p2, p3, flag)

    if nargin == 3   % Assume inhomogeneous coords
	flag = 'inhomog';
    end
    
    if ~all(size(p1)==size(p2)) | ~all(size(p1)==size(p3)) | ...
        ~(length(p1)==2 | length(p1)==3)                              
        error('points must have the same dimension of 2 or 3');
    end
    
    % If data is 2D, assume they are 2D inhomogeneous coords. Make them
    % homogeneous with scale 1.
    if length(p1) == 2    
        p1(3) = 1; p2(3) = 1; p3(3) = 1;
    end

    if flag(1) == 'h'
	% Apply test that allows for homogeneous coords with arbitrary
        % scale.  p1 X p2 generates a normal vector to plane defined by
        % origin, p1 and p2.  If the dot product of this normal with p3
        % is zero then p3 also lies in the plane, hence co-linear.
	r =  abs(dot(cross(p1, p2),p3)) < eps;
    else
	% Assume inhomogeneous coords, or homogeneous coords with equal
        % scale.
	r =  norm(cross(p2-p1, p3-p1)) < eps;
    end