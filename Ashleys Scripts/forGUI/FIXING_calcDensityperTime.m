function calcDensityperTime(obj,distances)
           %input:  distances in microns as an 1 x n vector of max values, values between are used
           %    example: distances = [100,200,300, inf]; interval is
           %    0<val<=100, 100<val<=200, 200<val<=300, val>300;
           % Average intensity within the distance is
           % make the distance mask
           
           if nargin<2
               distances = [5 40 200];
           end
           if isempty(obj.distmask)
               obj.makeDistanceMask();
           end
           if isempty(obj.cellFill.AIS_mask)
               regions = {'Total'};
           else
               regions = {'Total','No AIS','AIS only'};
           end
           %--- calculate background of SurfaceCargo Image
           max_cleanedcargomask = max(dip_image(obj.cleanedcargomask),[],3 %max project down to 2d imgae
           bgmaskthin = isnan(obj.distmask) & ~sum_cleanedcargomask; %get region without cleaned cargo mask and cellfill
           bgmask = berosion(bgmaskthin,5); %errode mask
           bgmask = repmat(bgmask,1,1,size(obj.cleanedcargomask,3)); %remake 3d
           
           fullbackgroundimage = GeneralAnalysis.regionfill_timeseries(obj.surfaceCargo.image*bgmask,~bgmask);% use regionfill to generate bg image
           SurfaceCargo_BGSubtracted = obj.surfaceCargo.image-fullbackgroundimage;
           SurfaceCargo_BGSubtracted(SurfaceCargo_BGSubtracted<0) = 0;
           
           currM.distance = distances;
           currM.rawintensity = zeros(numel(distances),size(SurfaceCargo_BGSubtracted,3));
           currM.areanormintensity = zeros(numel(distances),size(SurfaceCargo_BGSubtracted,3));
           for ii = 1:numel(distances)
               if ii == 1
                   currmask = obj.distmask<=distances(ii);
               else
                   currmask = obj.distmask>distances(ii-1) & obj.distmask<=distances(ii);
               end
               for rr = 1:numel(regions)
                   switch regions{rr}
                       case 'Total'
                           surfaceCargoMask = obj.cleanedcargomask;
                       case 'No AIS'
                           surfaceCargoMask = obj.cleanedcargomask-obj.cellFill.AIS_mask;
                       case 'AIS only'
                           surfaceCargoMask = obj.cleanedcargomask.*obj.cellFill.AIS_mask;
                   end
                   themask = currmask.*surfaceCargoMask;
                   signal = themask.*SurfaceCargo_BGSubtracted;
                   mask_area = sum(themask(:)); % area in number of pixels
                   % sliding window
                   windowsize = 4;
                   currM.rawintensity(ii,:) = movmean(signal,windowsize);
                   currM.areanormintensity(ii,:) = currM.rawintensity(ii,:)/mask_area; 
                   switch regions{rr}
                       case 'Total'
                           obj.M = currM;
                       case 'No AIS'
                           obj.M_noAIS = currM;
                       case 'AIS only'
                           obj.M_AIS = currM;
                   end
               end
      end
 end
              