classdef converter    
    methods(Static)
        function pupil_data = convert2method(pupil_data, method, data_mean, data_std)
            if strcmp(method, 'Z-score')
                pupil_data = converter.convert2z(pupil_data, data_mean, data_std);
                return;
            end
        end
        
        function z_data = convert2z(pupil_data, data_mean, data_std)
            z_data = (pupil_data-data_mean)/data_std;
        end
    end
end