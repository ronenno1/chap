classdef gui_lib
    methods(Static)
        function uic = uicontrol_checkbox(Figure_h, var, position)
            backgroundColor = [0.1 0.32 0.46];
            foregroundColor = [1 1 1];
            fontWeight      = 'bold';
            uic = uicontrol(Figure_h, 'Style', 'checkbox', ...
                                'BackgroundColor', backgroundColor, ...
                                'ForegroundColor', foregroundColor, ...
                                'Position', position,...
                                'FontWeight', fontWeight,...
                                'String', var);
        end

        function uisc = uicontrol_smart_checkbox(Figure_h, position, tooltipString, callback, keypressfcn)
            if (~exist('tooltipString', 'var'))
                tooltipString = '';
            end;
            if (~exist('callback', 'var'))
                callback = '';
            end;
            if (~exist('keypressfcn', 'var'))
                keypressfcn = '';
            end;
            foregroundColor = [1 1 1];
            backgroundColor = [0.1 0.32 0.46];

            uisc = uicontrol(Figure_h, 'Style','checkbox',...
                                      'BackgroundColor', backgroundColor,...
                                      'ForegroundColor', foregroundColor, ...
                                      'Position',position,...
                                      'Callback',callback,...
                                      'TooltipString', tooltipString,...
                                      'keypressfcn',keypressfcn);                       
        end
        
        function uit = uicontrol_title(string, x, y)
                uit = text(x, y, string,... 
                    'HorizontalAlignment', 'left',...
                    'Color', [0.07, 0.36, 0.55], ...
                    'FontWeight', 'bold',...
                    'FontSize',11);
        end

        function uit = uicontrol_text(Figure_h, string, position, tooltipString)
            if (~exist('tooltipString', 'var'))
                tooltipString = '';
            end;
            backgroundColor = [0.1 0.32 0.46];
            foregroundColor = [1 1 1];
            fontWeight      = 'bold';
            uit = uicontrol(Figure_h, 'Style', 'text', ...
                                'BackgroundColor', backgroundColor, ...
                                'ForegroundColor', foregroundColor, ...
                                'HorizontalAlignment', 'left',...
                                'Position', position,...
                                'String', string, ...
                                'FontWeight', fontWeight,...
                                'TooltipString', tooltipString);
        end
        
        function uie = uicontrol_edit(Figure_h, position, tooltipString, callback, keypressfcn)
            if (~exist('tooltipString', 'var'))
                tooltipString = '';
            end;
            if (~exist('callback', 'var'))
                callback = '';
            end;
            if (~exist('keypressfcn', 'var'))
                keypressfcn = '';
            end;
            uie = uicontrol(Figure_h, 'Style','Edit',...
                                      'BackgroundColor', [1 1 1],...
                                      'Position',position,...
                                      'Callback',callback,...
                                      'TooltipString', tooltipString,...
                                      'keypressfcn',keypressfcn);                       
        end
        
        function uib = uicontrol_button(Figure_h, position, string, callback, tooltipString)
            if (~exist('callback', 'var'))
                callback = '';
            end;
            if (~exist('tooltipString', 'var'))
                tooltipString = '';
            end;

            uib = uicontrol(Figure_h, 'Style', 'push',...
                                'String', string, ...
                                'Position', position, ...
                                'CallBack', callback,...
                                'ForegroundColor', [0.07, 0.36, 0.55], ...
                                'FontWeight', 'bold',...
                                'TooltipString', tooltipString);
            full_path = mfilename('fullpath');
            path      = full_path(1:strfind(full_path,mfilename)-1);

            imraw = imread([path 'images/buttun.jpg']);
            set(uib, 'unit', 'pixel');
            pos = get(uib, 'position');
            imfit = imresize(imraw, [pos(end) pos(end-1)]);
            set(uib, 'Cdata', imfit);            
        end
        
        function handler = paging(figure, handler, objs,  init_pos, max_obj2show, max_buts, prefix)
            if (~exist('prefix', 'var'))
                prefix = '';
            else
                prefix = strcat(prefix, '_');
            end;
            
            index = 0;
            init_pos_x = init_pos(1);
            init_pos_y = init_pos(2);
            w = init_pos(3);
            h = init_pos(4);
            text_pos_y = init_pos_y;
            for obj = objs
                index = index + 1;
                obj_text = [' ' strrep(char(obj), '_', ' ')];
                handler.(char(strcat(prefix, obj))) = gui_lib.uicontrol_checkbox(figure, obj_text, [init_pos_x text_pos_y w h]);
                text_pos_y = text_pos_y -h;
                if (index>max_obj2show)
                    set(handler.(char(strcat(prefix, obj))), 'Visible', 'off');
                end;
                if(mod(index, max_obj2show)==0)
                    text_pos_y = init_pos_y;
                end;
            end;
            if (index>max_obj2show)
                for i=1:min(max_buts, ceil(index/max_obj2show))
                    gui_lib.uicontrol_button(figure, [(i-1)*30+10 10 20 20], i, {@(~,~) gui_lib.go2(objs, handler, i, max_obj2show, prefix)});
                end;
            end;
        end
        
        function go2(header, handler, page, max_var2show, prefix)
            first = 1+max_var2show*(page-1);
            last  = max_var2show*page;
            for obj = 1:size(header, 2)
                visable = 'off';
                if(obj>=first && obj<= last)
                    visable = 'on';
                end;
                var = header(obj);
                set(handler.(char(strcat(prefix, var))), 'Visible', visable);
            end
        end

        
        function figure_h = create_figure(figure_name, background_path, dim)
            screensize = get(0,'ScreenSize');
            xpos = ceil((screensize(3)-dim(4))/2); % center the figure
            ypos = ceil((screensize(4)-dim(3))/2); % center the figure
            close(findobj('type','figure','name', figure_name))

            figure_h = figure('Name', figure_name,...
               'position',[xpos, ypos, dim(3), dim(4)],...
                'Color',[0.2 0.3 0.4],...
                'Resize','off',...
                'MenuBar','none',...
                'NumberTitle','off');
            gui_lib.set_background (background_path, dim);
        end

        function set_background (file_name, dim)
            full_path  = mfilename('fullpath');
            path       = full_path(1:strfind(full_path,mfilename)-1);
            background = imread([path file_name]);
            axes('Units','Pixels','Position',dim);
            image(background);
            axis off          % Remove axis ticks and numbers
            axis image        % Set aspect ratio to obtain square pixels
        end
        
        function val = get_checkbox_val(src)
            val  = get(src, 'Value');
        end

        function val = get_popupmenu_val(src)
            num  = get(src,'Value');
            vals = get(src,'String');
            val  = vals(num, :);
        end

    end
end

