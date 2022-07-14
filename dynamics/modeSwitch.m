function [newMode] = modeSwitch(mode, ie)

newMode = mode;
if (~isempty(ie))
    switch mode
        case 0 % dynamic grasping
            switch ie
                case 1
                    disp('Object losing contact')
                    newMode = -1;
                case 5
                    disp('Object slipping to the right')
                    newMode = 1;
                case 3
                    disp('Object slipping to the left')
                    newMode = 1;
                case 4
                    disp('Object tilting clockwise')
                    newMode = -1;
                case 2
                    disp('Object tilting counterclockwise')
                    newMode = -1;
                otherwise
            end
        case 1 % sliding
            switch ie
                case 1
                    disp('Object sliding off right side')
                    newMode = -1;
                case 2
                    disp('Object sliding off left side')
                    newMode = -1;
                case 3
                    disp('Object starting sticking again')
                    newMode = 0;
                otherwise
            end
        otherwise
    end
end



end

