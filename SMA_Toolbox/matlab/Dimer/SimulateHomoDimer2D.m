
classdef SimulateHomoDimer2D < handle
    % Dimer.SimulateHomoDimer2D -- This class simulates a 2D system of interacting
    % particles all of a single species A that undergo A+A <---> B reversible reactions
    % using a simmilar approach to Andrews & Bray particle-based stochastic simulation of
    % bimolecular reaction kinetics 
    % Ref: Andrews and Bray. Phys.Biol. 1(2004) 137--151
    %
    % Actually we are now basining the simulation off of the paper
    % Ref: Lipkova et. al., Analysis of Brownian Motion Dynamics Simulations of
    % Reverisible Bimolecular Reactions
    % This paper uses a reaction probability in addition to a reaction distance
    % to allow more accurate modeling of molecule position just before and after
    % dimerization events.

    % Concerns:
    %   1. The binding and unbinding radius timestep and diffusion constant
    %   must be tuned to get realistic kinetics
    %   2. This simulates the particle motion, but we should probably also
    %   simulate the particle observation on a coarser timescale.  The
    %   simulation is more accurate when timesteps are small, but if the
    %   observation and simulation timescales are synced the more accurate we
    %   make the simulation, the timescales may no longer be realistic for
    %   fluorecence microscopy which has practiocal finite time resolution.
    %   3.  The Andrews&Bray algorithm uses a unbinding step where the particles are
    %   forced appart to a given fixed radius with a sudden increase in
    %   displacement.  This part of the simulation is only accurate at large
    %   space and time scales, and may lead to problems with dimerization
    %   detection.  e.g. it may present artifical signals to the estimation
    %   algorithm which can be picked up to detect the dissociation of a dimer.
    %   To avoid this we need to make the simulation time steps small, which
    %   allows the unbinding radius to be smaller.
    
    properties
        %These parameters are set in the constructor
        nparticles; %Number of particles to simulate
        nsteps; %number of total simulation steps to make

        %These parameters can be modified before simulation to control the model of the reactions
        simulation_subsampling=8;  % number of sub-steps to take for each simulated observation (at least 100 for accurate motion blur)
        box_size=10.0; % (microns) the size of the box the particles are constrained to
        binding_radius = 0.050 %(microns) Distance at which a dimer will form between 2 interacting particles
        binding_probability = 0.20 % (1/frame) Probility of binding per second when close enough.
        unbinding_radius = 0.150 %(micron) Distance to which individual paricles will be forced when a dimer is broken.
        koff = 5.0e-1; % Dissociation rate (1/s)
        monomer_diffusion_constant = 0.050% um^2/s
        dimer_diffusion_constant   = 0.020% um^2/s

        %The position sigmas will be modeled as gamma distributed
        pos_sigma_mean = 0.030; % microns
        pos_sigma_stddev = 0.005; %microns
        pos_sigma_shape = 6; % This is the shape parameter for the gamma distribution 
        frame_time = 0.05; % seconds

        %computed parameters
        sim_time_step; % (seconds) (computed from frame_time and simulation_subsampling
        max_time; %time to simulate for this is compuited from sim_time_
        
        tracks; %The particle tracks: size:[nparticles, nsteps, ndim]
        times;  %The times for each step: size:[1,nsteps]
        dimer_sequence; 
        associations; % size:[nassociations,3] cols are timestep, pair_i, pair_j
        nassociations=0;
        dissociations; % size:[nassociations,3] cols are timestep, pair_i, pair_j
        ndissociations=0;
    end

    methods 
        function obj = SimulateHomoDimer2D(nparticles_, nsteps_)
            obj.nparticles = nparticles_;
            obj.nsteps = nsteps_;
            obj.sim_time_step = obj.frame_time/obj.simulation_subsampling;
            obj.max_time = (obj.nsteps-1) * obj.sim_time_step;          
        end
        

%         function pairid=findDimer(obj)
%         end
% 
%         function [pair_mat, assocs, dissocs, state] = makePair(obj,pairid)
%             if pairid(1)>pairid(2)
%                 pairid = pairid(2:1);
%             end
%             assert(pairid(1)~=pairid(2));
%             pair_mat=zeros(PairAnalysis.NPairColumns, obj.n
%         end
        
%         function [pair_id, pair_mat, assocs, dissocs, state] = getInteractingPair(obj)
%             nassoc_events = zerso(obj.nparticles, obj.nparticles);
%             for k=1:nassociations
%                 [i,j] = associations(i,[2,3]);
%                 nassoc_events(i,j)=nassoc_events(i,j)+1;
%             end
%             [max_nassoc_events, max_pair_idx] = max(nassoc_events(:));
%             if max_nassoc_events==0
%                 error('SimulateHomoDimer2D:getInteractingPair','No associations detected');
%             end
%             pair_id = ind2sub([ obj.nparticles, obj.nparticles],  max_pair_idx);
%             
%             tracks = obj.tracks(pair,:,:);
%             
%             
%         end
        
        function pairids = identifyMostReactivePairs(obj, max_pairs)
            if nargin==1
                max_pairs=0;
            end
            interactions = zeros(obj.nparticles);
            for i=1:size(obj.associations,1)
                interactions(obj.associations(i,2), obj.associations(i,3)) =  interactions(obj.associations(i,2), obj.associations(i,3))+1;                
            end
            [vals,sidx] = sort(interactions(:),'descend');
            n_int_pairs = find(vals==0,1,'first')-1; %number of interacting pairs
            if max_pairs>0
                n_int_pairs = min(n_int_pairs,max_pairs);
            end
            pairids = zeros(n_int_pairs, 2);
            for i=1:n_int_pairs
                [n,m]=ind2sub(size(interactions),sidx(i));
                pairids(i,:) = [n,m];
            end                               
        end

        function [pair_mat, assocs, dissocs, state, pairid]=observePair(obj,pairid)
            if nargin==1
                pairid = obj.identifyMostReactivePairs(1);
            end
            assert(numel(pairid)==2)
            if pairid(1)>pairid(2)
                pairid = pairid([2,1]);
            end
            nsamples = floor(obj.nsteps/obj.simulation_subsampling);
            sample_times = (0:nsamples-1)'*obj.frame_time;
            
            pair_mat=zeros(nsamples,PairAnalysis.NPairColumns);
            pair_mat(:,1) = sample_times;
            %Gamma distributed SE's for positions of both particles
%             pair_mat(:,[10,11,14,15]) = gamrnd(obj.pos_sigma_shape, obj.pos_sigma_mean/obj.pos_sigma_shape, nsamples, 4);
            %Sample as normally distributed
            pair_mat(:,[10,11,14,15]) = abs(obj.pos_sigma_mean+obj.pos_sigma_stddev*randn(nsamples, 4));
            pair_mat(:,18) = (1:nsamples)'; %sample times
            
            %Get motion blurred mean positions
            for si = 1:nsamples
                idxs = (((si-1)*obj.simulation_subsampling):(si*obj.simulation_subsampling-1))+1;
                %particle A
                pair_mat(si,2) = mean(squeeze(obj.tracks(pairid(1),idxs,1)));
                pair_mat(si,3) = mean(squeeze(obj.tracks(pairid(1),idxs,2)));
                %particle B
                pair_mat(si,6) = mean(squeeze(obj.tracks(pairid(2),idxs,1)));
                pair_mat(si,7) = mean(squeeze(obj.tracks(pairid(2),idxs,2)));
            end
            %simulate measurment errors
            pair_mat(:,2) = normrnd(pair_mat(:,2),pair_mat(:,10));
            pair_mat(:,3) = normrnd(pair_mat(:,3),pair_mat(:,11));
            pair_mat(:,6) = normrnd(pair_mat(:,6),pair_mat(:,14));
            pair_mat(:,7) = normrnd(pair_mat(:,7),pair_mat(:,15));            
            
            if obj.nassociations==0
                assocs=[];
            else
                assocs = obj.associations( obj.associations(:,2)==pairid(1) & obj.associations(:,3)==pairid(2), 1);
            end
            
            if obj.ndissociations==0
                dissocs=[];
            else
                dissocs = obj.dissociations( obj.dissociations(:,2)==pairid(1) & obj.dissociations(:,3)==pairid(2), 1);
            end
            changepoints = sort([assocs;dissocs]);
            
            state=false(1,nsamples);
            for i=1:numel(changepoints)
                idx = ceil(changepoints(i)/obj.simulation_subsampling);
                state(idx:end) = ~state(idx:end); 
            end
            assocs = assocs * obj.sim_time_step;
            dissocs = dissocs * obj.sim_time_step;
        end
        
        function plotPairMatrixDistances(obj, pairid)
            if nargin==1
                pairid = obj.identifyMostReactivePairs(1);
            end
            [pair, assocs, dissocs, state]=obj.observePair(pairid);
            f=figure();
            certainty=0.95;
            D = PairAnalysis.makePairDists(pair, certainty);
            plot(D(:,1),D(:,2),'-r','DisplayName','Dist');
            hold on;
            true_dists = sqrt(sum((squeeze(obj.tracks(pairid(1),:,:)) - squeeze(obj.tracks(pairid(2),:,:))).^2,2));
            plot( (0:obj.nsteps-1)'*obj.sim_time_step,true_dists,':r','LineWidth',0.5,'DisplayName','True Dist.');
%             plot(D(:,1),D(:,3),'--r','DisplayName',sprintf('Min Dist (%.1f%% conf.)',certainty*100));
%             plot(D(:,1),D(:,4),'--r','DisplayName',sprintf('Max Dist (%.1f%% conf.)',certainty*100));
            xs = [D(:,1)', fliplr(D(:,1)')];
            ys = [D(:,3)', fliplr(D(:,4)')];
            fill(xs,ys,'r','EdgeColor','none','FaceAlpha',0.4,'DisplayName',sprintf('Dist (%.1f%% conf.)',certainty*100));
            yl = ylim();
%             for i=1:numel(state)
%                 xs=[D(i,1),D(i,1)+obj.frame_time, D(i,1)+obj.frame_time, D(i,1)];
%                 ys=[0,0, yl(2), yl(2)];
%                 if state(i)
%                     h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
%                     h.Annotation.LegendInformation.IconDisplayStyle='off';
%                 end
%                 xs=[D(i,1),D(i,1)];
%                 ys=[0, yl(2)];
%                 h=fill(xs,ys,'-','EdgeColor','k', 'EdgeAlpha',0.3, 'FaceColor','none');
%                 h.Annotation.LegendInformation.IconDisplayStyle='off';
%             end
            
            plot([D(1,1),D(end,1)], [obj.binding_radius obj.binding_radius], ':k','DisplayName', 'Binding Radius');
            plot([D(1,1),D(end,1)], [obj.unbinding_radius obj.unbinding_radius], '--k','DisplayName','Unbinding Radius');
            for i=1:numel(assocs)
                h=plot([assocs(i), assocs(i)],[0,yl(2)],'c-','LineWidth',2);
                h.Annotation.LegendInformation.IconDisplayStyle='off';
            end
            for i=1:numel(dissocs)
                h=plot([dissocs(i), dissocs(i)],[0,yl(2)],'m-','LineWidth',2);
                h.Annotation.LegendInformation.IconDisplayStyle='off';
            end
            for i=1:numel(state)
            end
            hold off;
            ylim(yl);
            xlabel('Time (s)');
            ylabel('Dist (um)');
            legend('Location','best');
            
        end
        
        function plotPairObservations(obj, pairid, show_simulation)
            if nargin==1
                pairid = obj.identifyMostReactivePairs(1);
            end
            if nargin<=2
                show_simulation=1;
            end
            [pair_mat, assocs, dissocs, state]=obj.observePair(pairid);
            f=figure();
            p1=min(pairid);
            p2=max(pairid);
            ax=gca();
            hold on;
            markersize=3;
            plot3(pair_mat(:,2), pair_mat(:,3), pair_mat(:,1),'or-','MarkerSize',markersize,'LineWidth',2.0,'MarkerFaceColor',[1,0.3,0.3]);
            plot3(pair_mat(:,6), pair_mat(:,7), pair_mat(:,1),'ob-','MarkerSize',markersize,'LineWidth',2.0,'MarkerFaceColor',[0.3,0.3,1]);
            alpha(0.5);
            if show_simulation
                plot3(obj.tracks(p1,:,1), obj.tracks(p1,:,2), obj.times,'-', 'LineWidth',1.0,'color',[0.5,0,0]);
                plot3(obj.tracks(p2,:,1), obj.tracks(p2,:,2), obj.times,'-', 'LineWidth',1.0,'color',[0,0,0.5]);
                for i=1:length(assocs)
                    at=assocs(i)/obj.sim_time_step;
                    if i>length(dissocs)
                        dt=obj.nsteps;
                    else
                        dt=dissocs(i)/obj.sim_time_step;
                    end
                    idxs = (at:dt);
                    xs=0.5*(squeeze(obj.tracks(p1,idxs,1)+obj.tracks(p2,idxs,1)));
                    ys=0.5*(squeeze(obj.tracks(p1,idxs,2)+obj.tracks(p2,idxs,2)));
                    plot3(xs,ys,obj.times(idxs),'-','color',[0.0,1.0,0.5],'linewidth',3.0);                
                end
            end
            hold off;
            ax.Box='on';
            ax.BoxStyle='Full';
            ax.XGrid='on';
            ax.YGrid='on';
            ax.ZGrid='on';
            ax.XMinorGrid='on';
            ax.YMinorGrid='on';
            ax.ZMinorGrid='on';
            pixel_size=0.1; %micron
%             ax.XTick=(0:obj.box_size/pixel_size)*pixel_size;
%             ax.YTick=(0:obj.box_size/pixel_size)*pixel_size;            
            ax.GridAlpha=0.5;
            xl=xlim();
            yl=ylim();
            yratio = (xl(2)-xl(1))/(yl(2)-yl(1));
            pbaspect([yratio 1 2])
            xlabel('X (um)');
            ylabel('Y (um)');
            zlabel('Time (s)');
        end
% 
        function plotParticles(obj)
            f=figure(1);
            set(f,'Renderer','OpenGL')
            clf;
            colormap('jet') 
            cmap=colormap;
            cmapsize=int64(64);
            ax=axes();
            hold on;
            for n=1:obj.nparticles
                idx=(mod(n,min(cmapsize,obj.nparticles))+1) * cmapsize/(min(cmapsize,obj.nparticles));
                plot3(obj.tracks(n,:,1), obj.tracks(n,:,2), obj.times,'Color',cmap(idx,:));
            end
            hold off;
            ax.Box='on';
            ax.BoxStyle='Full';
            ax.XGrid='on';
            ax.YGrid='on';
            ax.ZGrid='on';
            ax.XMinorGrid='on';
            ax.YMinorGrid='on';
            ax.ZMinorGrid='on';
            xl=xlim();
            yl=ylim();
            yratio = (xl(2)-xl(1))/(yl(2)-yl(1));
            pbaspect([1 1/yratio 1])
            xlabel('X (um)')
            ylabel('Y (um)')
            zlabel('Time (s)')
        end

        function plotAnimatedParticles(obj)
            f=figure(1);
            set(f,'Renderer','OpenGL')
            clf;
            colors = jet(obj.nparticles);
            sizes = repmat(30,1,obj.nparticles);
            pbaspect([1 1 1])
            for t=1:obj.nsteps
                clf;
                axis([0 obj.box_size 0 obj.box_size]);
                hold on;
                for di=find(obj.dimer_sequence(:,t)') % Force dimers to have same dx movement
                    partner=obj.dimer_sequence(di,t);
                    if di<partner %Only the smaller of the two indicies in a pair will call this
                        line([obj.tracks(di,t,1) obj.tracks(partner,t,1)], ...
                             [obj.tracks(di,t,2) obj.tracks(partner,t,2)], 'Color','r','LineWidth', 4.0); 
                    end
                end
                scatter(obj.tracks(:,t,1),obj.tracks(:,t,2),sizes',colors,'filled');
                hold off;
                title(sprintf('Frame:%i',t));
                drawnow; pause(0.001);
            end
           
        end


        function simulateTracks(obj)
            dimers=zeros(1,obj.nparticles); % if (i,j) is a dimer and i<j dimer(i)=j, dimer(j)=i, otherwise dimer(i)=0.
            obj.tracks=zeros(obj.nparticles,obj.nsteps,2); % dims are: particle# x timestep# x [x,y] pos
            obj.associations=zeros(obj.nparticles,3); %preallocate
            obj.dissociations=zeros(obj.nparticles,3); %preallocate
            obj.nassociations=0;
            obj.ndissociations=0;
            obj.dimer_sequence=zeros(obj.nparticles, obj.nsteps);
            
            obj.tracks(:,1,:)=rand(obj.nparticles,2)*obj.box_size;  %starting positions
            for step=2:obj.nsteps
                diffusion_constant=obj.monomer_diffusion_constant.*(dimers==0)' + ...
                                    obj.dimer_diffusion_constant.*(dimers~=0)';
                diffusion_constant=repmat(diffusion_constant,1,2);
                dx=sqrt(2*obj.sim_time_step*diffusion_constant).*randn(obj.nparticles,2); %displacements 
                for di=find(dimers) % Force dimers to have same dx movement
                    if di<dimers(di) %Only the smaller of the two indicies in a pair will call this
                        
                        dx(di,:)=dx(dimers(di),:);
                    end
                end
                positions=squeeze(obj.tracks(:,step-1,:))+dx;
                for n=1:obj.nparticles  %reflect at the boundaries
                    if positions(n,1)<0
                        positions(n,1)=-positions(n,1);
                    end
                    if positions(n,1)>obj.box_size
                        positions(n,1)=obj.box_size-(positions(n,1)-obj.box_size);
                    end
                    if positions(n,2)<0
                        positions(n,2)=-positions(n,2);
                    end
                    if positions(n,2)>obj.box_size
                        positions(n,2)=obj.box_size-(positions(n,2)-obj.box_size);
                    end
                end
                for di=find(dimers) % Force dimers to have same dx movement
                    if di<dimers(di) %Only the smaller of the two indicies in a pair will call this
                        partner=dimers(di);
                        [np1, np2]=obj.setDimerPositions(positions(di,:), positions(partner,:));
                        positions(di,:)=np1;
                        positions(partner,:)=np2;
                    end
                end

               
                old_dimers=dimers;
                %Make dimers
                for n=1:obj.nparticles
                    if dimers(n) %This particle is already in a dimer
                        continue
                    end
                    for m=n+1:obj.nparticles
                        if dimers(m) %This particle is already in a dimer
                            continue
                        end
                        vec=positions(m,:)-positions(n,:);
                        d=sqrt(sum(vec.*vec)); %distance
                        if ( d < obj.binding_radius ) && ( rand < obj.binding_probability )
                            [np1, np2]=obj.setDimerPositions(positions(n,:), positions(m,:));
                            positions(n,:)=np1;
                            positions(m,:)=np2;
                            dimers(m)=n;
                            dimers(n)=m;
                            obj.recordAssociation(step,n,m);
                        end
                    end
                end
                %Break dimers
                for di=find(old_dimers)
                    if di<old_dimers(di) %Only the smaller of the two indicies in a pair will call this
                        if -log(rand ) < obj.sim_time_step*obj.koff;
                            phi=2.*pi*rand;
                            v=[cos(phi) sin(phi)];
                            partner=old_dimers(di);
                            cmass=positions(di,:);
                            positions(di,:)=cmass+obj.unbinding_radius*0.5*v;
                            positions(partner,:)=cmass-obj.unbinding_radius*0.5*v;
                            dimers(di)=0;
                            dimers(partner)=0;
                            obj.recordDissociation(step,di,partner);
                        end
                    end
                end
                %update final positions
                obj.tracks(:,step,:)=positions;                
                obj.dimer_sequence(:,step)=dimers;
            end
            obj.times=((1:obj.nsteps)-1)*obj.sim_time_step;
            if obj.nassociations==0
                obj.associations=[];
            else
                obj.associations(obj.nassociations+1:end,:)=[];
            end
            if obj.ndissociations==0
                obj.dissociations=[];
            else
                obj.dissociations(obj.ndissociations+1:end,:)=[];
            end
        end
    end
    methods (Access=private)
        function recordAssociation(obj,step,n,m)
            s=length(obj.associations);
            if obj.nassociations==s %double size
                new_assocs=zeros(s*2,3);
                new_assocs(1:s,:)=obj.associations;
                obj.associations=new_assocs;
            end
            obj.nassociations=obj.nassociations+1;
            obj.associations(obj.nassociations,:)=[step n m];
        end
        
        function recordDissociation(obj,step,n,m)
            s=length(obj.dissociations);
            
            if obj.ndissociations==s %double size
                new_dissocs=zeros(s*2,3);
                new_dissocs(1:s,:)=obj.dissociations;
                obj.dissociations=new_dissocs;
            end
            obj.ndissociations=obj.ndissociations+1;
            obj.dissociations(obj.ndissociations,:)=[step n m];
        end

        function [np1,np2]=setDimerPositions(obj, p1, p2)
            cmass=0.5*(p1+p2);
            phi=2.*pi*rand;
            vec=[cos(phi) sin(phi)];
            np1=cmass-vec*0.5*obj.binding_radius;
            np2=cmass+vec*0.5*obj.binding_radius;
            while ~obj.inBounds(np1) || ~obj.inBounds(np1)
                if ~obj.inBounds(np1)
                    [np1, np2]=obj.setDimerPositions(np2,np2);
                elseif ~obj.inBounds(np2)
                    [np1, np2]=obj.setDimerPositions(np1,np1);
                end
            end
        end

        function ok=inBounds(obj,p)
            if p(1)<0 || p(1)>obj.box_size || p(2)<0 || p(2)>obj.box_size
                ok=false;
            else
                ok=true;
            end
        end
    end
end
