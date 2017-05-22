classdef MTCell <handle
    %% After creating MT cell, before using with V1 Cells, according to direction selectivity of V1, 
    % use MTCell.setKernel(V1Selectivity)
    properties(Constant)
        NO_SURROUND = 'noSurround';
        SYMMETRIC_ISOTROPIC = 'symmetricIsotropic';
        SYMMETRIC_ANISOTROPIC = 'symmetricAnisotropic';
        ASYMMETRIC_ANISOTROPIC = 'asymmetricAnisotropic';
    end
    
    properties(Access = 'public')
        E_exc = 70;
        E_inh = -10;
        E_l = 0;
        orientation,
        density, 
        receptiveFieldSize,
        kernel,
        responseFrames
        
        CreateMethod
        
    end
    
    methods
        function this = MTCell(orientation, d, CreateMethod, receptiveFieldSize, leak)
            
            if ~exist('orientation', 'var') || isempty(orientation)
                orientation = 0;
            end
            
            if ~exist('d', 'var') || isempty(d)
                d = 0.3;
            end
            
            if ~exist('receptiveFieldSize', 'var') || isempty(receptiveFieldSize)
                receptiveFieldSize = 20;
            end
            
            if ~exist('leak', 'var') || isempty(leak)
                leak = 0.1;
            end
    
            this.orientation = orientation;
            this.density = d;
            this.receptiveFieldSize = receptiveFieldSize;
            this.E_l = leak;
            this.setCreateMethod(CreateMethod);
        end
    end
    
    methods(Access = 'public')
        function out = Integrate(this, V1ResponsesAlongFrames, duration)
            this.responseFrames = [];
            for i = 1:duration:size(V1ResponsesAlongFrames,3)
                responseFrame = this.IntegrateAlongDuration(V1ResponsesAlongFrames, i, min(i+duration-1, size(V1ResponsesAlongFrames,3)));
                minV = min(min(responseFrame));
                maxV = max(max(responseFrame));
                if(maxV - minV>0)
                    responseFrame = (responseFrame - minV)./(maxV - minV);
                end
                this.responseFrames = cat(3, this.responseFrames, responseFrame);
            end
            out = this.responseFrames;
        end
        
        function out = IntegrateAlongDuration(this, V1ResponsesAlongFrames, start, ending)
            rFrames = [];
            for i = start:ending
                rFrame = calculateResponse(this, V1ResponsesAlongFrames);
                rFrames = cat(3, rFrames, rFrame);
            end
            
            out = sum(rFrames, 3);
        end
        
        function out = calculateResponse(this, V1ResponsesAlongFrames)
            coeffs = zeros(size(this.kernel));
            center = ceil(size(this.kernel)./2);
            coeffs(:,:) = this.E_inh;
            coeffs(center(1)-floor(this.receptiveFieldSize/4): center(1)+floor(this.receptiveFieldSize/4), center(2)-floor(this.receptiveFieldSize/4):center(2)+floor(this.receptiveFieldSize/4)) = this.E_exc;
            temp_kernel = this.kernel.*coeffs;
            if size(V1ResponsesAlongFrames,3)>0
                kernel3d = repmat(temp_kernel./size(V1ResponsesAlongFrames,3), 1, 1, size(V1ResponsesAlongFrames,3));
            else
                kernel3d = repmat(temp_kernel, 1, 1, size(V1ResponsesAlongFrames,3));
            end
            response = imfilter(V1ResponsesAlongFrames, kernel3d, 'conv');
            response = response./(this.E_exc + this.E_inh + this.E_l );
            response = sum(response,3);
            out = response;
        end
        
        function setKernel(this, v1_orient)
            radius = this.receptiveFieldSize;
            distances_x = repmat(-floor(radius/2):1:floor(radius/2), radius,1);
            distances_y = repmat((-floor(radius/2):1:floor(radius/2)), radius,1);
            distances = sqrt(distances_x.^2 + distances_y.^2);
            clearvars distances_x distances_y
            sigma_c = this.receptiveFieldSize/3;
            sigma_s = 4*this.receptiveFieldSize/3;
            sigma_ss = this.receptiveFieldSize/3;
            wc = exp(-distances.^2./(2*sigma_c^2));%/(sigma_c*sqrt(2*pi));
            wc = wc./(sigma_c*sqrt(2*pi));
            
            if strcmpi(this.CreateMethod, MTCell.NO_SURROUND)
                ws = 0;
            elseif strcmpi(this.CreateMethod, MTCell.SYMMETRIC_ISOTROPIC)
                ws = exp(-distances.^2./(2*sigma_s^2));
                ws = ws./(sigma_s*sqrt(2*pi));
            elseif strcmpi(this.CreateMethod, MTCell.SYMMETRIC_ANISOTROPIC)
                angle_mt = this.orientation;
                u = [radius/2 + 3*sigma_ss*cos(angle_mt); radius/2 + 3*sigma_ss*sin(angle_mt)];
                v = [radius/2 + 3*sigma_ss*cos(angle_mt+pi); radius/2 + 3*sigma_ss*sin(angle_mt+pi)];
                temp_u = (distances-u(1)).*(distances-u(1)) + (distances-u(2)).*(distances-u(2));
                temp_v = (distances-v(1)).*(distances-v(1)) + (distances-v(2)).*(distances-v(2));
                ws = exp(-temp_u./(2*sigma_ss^2)) + exp(-temp_v./(2*sigma_ss^2));
                ws = ws./(sigma_ss*sqrt(2*pi));
            elseif strcmpi(this.CreateMethod, MTCell.ASYMMETRIC_ANISOTROPIC)
                angle_mt = this.orientation;
                u = [radius/2 + 3*sigma_ss*cos(angle_mt); radius/2 + 3*sigma_ss*sin(angle_mt)];
                temp_u = (distances-u(1)).*(distances-u(1)) + (distances-u(2)).*(distances-u(2));
                ws = exp(-temp_u./(2*sigma_ss^2));
                ws = ws./(sigma_ss*sqrt(2*pi));
            end
            wd = wc-ws;
            if abs(this.orientation - v1_orient)<pi/2
                this.kernel = cos(abs(this.orientation - v1_orient))*wd;
            elseif abs(this.orientation - v1_orient)>pi/2
                this.kernel = -1*cos(abs(this.orientation - v1_orient))*wd;
            else
                this.kernel = wd*0;
            end
            

        end
    end
    
    methods(Access = 'private')
        function setCreateMethod(this, createMethod)
            if ~exist('createMethod', 'var') || isempty(createMethod)
                this.CreateMethod = MTCell.NO_SURROUND;
            elseif strcmpi((createMethod), (MTCell.NO_SURROUND))
                this.CreateMethod = MTCell.NO_SURROUND;
            elseif strcmpi((createMethod), (MTCell.SYMMETRIC_ISOTROPIC))
                this.CreateMethod = MTCell.SYMMETRIC_ISOTROPIC;
            elseif strcmpi((createMethod), (MTCell.SYMMETRIC_ANISOTROPIC))
                this.CreateMethod = MTCell.SYMMETRIC_ANISOTROPIC;
            elseif strcmpi((createMethod), (MTCell.ASYMMETRIC_ANISOTROPIC))
                this.CreateMethod = MTCell.ASYMMETRIC_ANISOTROPIC;
            end
            
        end
        
    end
end