classdef LDhs < geometry.earthModel
    properties
        % rigidity
        G;
        % Poisson's ratio
        nu;
    end
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function o=LDhs(G,nu)
            % LDhs is a class providing the necessary functions to
            % compute the stress interactions between sources and receivers
            % and among the receivers.
            %
            %   earthModel = LDhs(G,nu);
            %
            % where G is rigidity and nu is the Poisson's ratio of a
            % homogeneous elastic half space for 2-d line displacement.
            %
            
            if (0==nargin)
                return
            end
            
            assert(0<=G,'LDhs::rigidity must be positive.')
            assert(nu<=0.5,'LDhs::Poisson''s ratio should be lower than 0.5.')
            assert(-1<=nu,'LDhs::Poisson''s ratio should be greater than -1.')
            
            o.G=G;
            o.nu=nu;
        end
        
        function [varargout]=tractionKernels(obj,src,rcv)
            % TRACTIONKERNELS computes the traction on receiver faults due
            % to motion of rectangular dislocations in a half space.
            %
            % rcv - receiver fault
            %
            % SEE ALSO: unicycle
            
            
            varargout = cell(1,nargout);
            [varargout{:}]=computeTractionKernels(src,rcv);
        end
        
        function [varargout]=stressKernels(obj,rcv)
            % STRESSKERNELS computes the stress on receiver shear zone due
            % to motion of rectangular dislocations in a half space.
            %
            % rcv - receiver fault
            %
                        
            
            varargout = cell(1,nargout);
            [varargout{:}]=computeStressKernels(obj,rcv);
        end
        
        function [varargout]=displacementKernels(obj,src,obs)
            % DISPLACEMENTKERNELS computes the stress on receiver faults due to
            % motion of rectangular dislocations in a half space.
            %
            % src - source fault
            %
            
            
            varargout = cell(1,nargout);
            [varargout{:}]=computeDisplacementKernels(src,obs);
        end
    end % methods
end % class definition