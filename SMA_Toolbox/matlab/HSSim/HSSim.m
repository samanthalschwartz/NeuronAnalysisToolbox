classdef HSSim
    %HSSIM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        imageSize;
        spectralBounds=[400 750];
        psfSigma=[1, 1.14];
        spotSize2D;
        spotSize3D;
        spotIm2D;
        spotIm3D;
        meanI=[500 1000 1500 2000]; %Different I means
        meanIdist=[0.4 0.4 0.1 0.1]; %Proabability of selecting a particular meanI
        meanBG=3;
        stdBG=1;
    end
    
    properties (Dependent=true)
        imageSize2D
        imageSize3D
    end

    methods
        function obj=HSSim(imageSize, psfSigma)
            obj.imageSize=imageSize;
            obj.psfSigma=psfSigma;
            if length(obj.imageSize)==2
                obj.imageSize(3)=1;
            end
            if length(obj.psfSigma)==2
                obj.psfSigma(3)=1;
            end

            spotSize=ceil(2*3*obj.psfSigma+1);
            spotSize=spotSize+mod(spotSize+1,2); %Make it odd
            obj.spotSize2D=spotSize(1:2);
            obj.spotSize2D=obj.spotSize2D(:)';
            obj.spotSize3D=spotSize(1:3);
            obj.spotSize3D=obj.spotSize3D(:)';

            obj.spotIm2D=gaussf(deltaim(obj.spotSize2D),obj.psfSigma(1:2));
            if obj.imageSize(3)>1
                obj.spotIm3D=newim(obj.spotSize3D);
                for n=1:obj.spotSize3D(3)
                    norm=1/(sqrt(2*pi)*obj.psfSigma(3))*exp(-(n-ceil(obj.spotSize3D(3)/2))^2/(2*obj.psfSigma(3)));
                    im=gaussf(norm*deltaim(obj.spotSize2D),obj.psfSigma(1:2));
                    obj.spotIm3D(:,:,n-1)=im;
                end
            end
        end

        function intensity=randIntensity(obj)
            if isscalar(obj.meanI)
                m=meanI;
            else
                r=rand();
                for i=1:length(obj.meanI)
                    if r<=obj.meanIdist(i)
                        %select this index as mean
                        m=obj.meanI(i);
                        break;
                    end
                    r=r-obj.meanIdist(i);
                end
            end
            stdI=sqrt(m); %approximate std deviation of particle intensity as a function of the mean
            intensity=randn(1,1)*stdI+m;
        end
        
        function [outim,outpos]=simulate2DImage(obj,nParticles,nT)
            if nargin<3
                nT=1;
            end
            spotDisp=floor(obj.spotSize2D/2);
            outim=newim([obj.imageSize2D nT]);
            outpos=zeros(nParticles*nT,3);
            for t=1:nT
                pos=floor(rand(nParticles,2).*repmat(obj.imageSize2D,nParticles,1));
                im=newim(obj.imageSize2D+2*spotDisp);
                for n=1:nParticles
                    tl=pos(n,:);
                    br=pos(n,:)+obj.spotSize2D-1;
                    intensity=obj.randIntensity();
                    im(tl(1):br(1), tl(2):br(2))=im(tl(1):br(1), tl(2):br(2))+intensity*obj.spotIm2D;
                end
                outim(:,:,t-1)=im(spotDisp(1):end-spotDisp(1),spotDisp(2):end-spotDisp(2));
                outpos(nParticles*(t-1)+1:nParticles*t,1:2)=pos+1;
                outpos(nParticles*(t-1)+1:nParticles*t,3)=t;
            end
            outim=single(obj.addNoise(outim));
        end

        function [outim,outpos]=simulate3DImage(obj,nParticles,nT)
            if nargin<3
                nT=1;
            end
            spotDisp=floor(obj.spotSize3D/2);
            outim=newim([obj.imageSize3D nT]);
            outpos=zeros(nParticles*nT,4);
            for t=1:nT
                pos=floor(rand(nParticles,3).*repmat(obj.imageSize3D,nParticles,1));
                im=newim(obj.imageSize3D+2*spotDisp);
                for n=1:nParticles
                    tl=pos(n,:); %top left of image tile
                    br=pos(n,:)+obj.spotSize3D-1; %bottom right of image tile
                    intensity=obj.randIntensity();
                    im(tl(1):br(1), tl(2):br(2), tl(3):br(3))=im(tl(1):br(1), tl(2):br(2), tl(3):br(3))+intensity*obj.spotIm3D;
                end
                outim(:,:,:,t-1)=im(spotDisp(1):end-spotDisp(1),spotDisp(2):end-spotDisp(2),spotDisp(3):end-spotDisp(3));
                outpos(nParticles*(t-1)+1:nParticles*t,1:3)=pos;
                outpos(nParticles*(t-1)+1:nParticles*t,4)=t;
            end
            outim=single(obj.addNoise(outim));
        end

        function val=get.imageSize2D(obj)
            val=obj.imageSize(1:2);
        end
        function val=get.imageSize3D(obj)
            val=obj.imageSize(1:3);
        end

        function im=addNoise(obj, im)
            nim=abs(noise(newim(size(im))+obj.meanBG,'gaussian',obj.stdBG));
            im=abs(round(im+nim));
            im=noise(im,'poisson',1);
        end

    end %Public methods 
end %Classdef

