% %%Clear variables
% function cleariface()
%     if exist('iface')
%         try
%             Gauss2DModel_iface('@delete', iface);
%         catch err
%             if strcmp(err.identifier, 'ClassHandle:covertMat2HandlePtr')
%                 clear('iface'); %Already deleted but we still have handle
%             end
%         end
%     end
% 
% 
%     mex_name='Gauss2DModel_iface';
%     mex_module=getIfaceMexModule(mex_name);
%     if ~isempty(mex_module)
%         fprintf('Clearing: %s\n',mex_module);
%         clear(mex_module);
%     else
%         fprintf('No module "%s" found\n', mex_name);
%     end
% 
% end
% 
% function mex_module=getIfaceMexModule(mex_name)
%     mex_name=[mex_name '.' mexext];
%     [~,f] = inmem( '-completenames' );
%     l=not(cellfun(@isempty,strfind(f,mex_name)));
%     if any(l)
%         mex_module= f{ l };
%     else
%         mex_module='';
%     end
% end
