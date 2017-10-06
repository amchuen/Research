classdef fieldVar < handle
% Abstract class for creating field variables when solving PDEs
% Subclass constructor should accept data for the field values,
% boundary conditions, and the control switches for iterative methods.	
% Should be able to implement vector calculus through numerical methods once field variable is
% implemented as a subclass.
    
    properties (Abstract)
        % setup
        gr % grid
        fv % field value
        bc % boundary conditions
        ct % control variables
       
    end
    
    methods (Abstract)
        % All field variables should be able to do the following:
        % Take the laplacian
        % Take the gradient/jacobian
        % Take the Divergence
        
        
    end
    
end