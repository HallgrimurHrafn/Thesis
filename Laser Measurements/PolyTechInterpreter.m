function Output=PolyTechInterpreter(freq,volt)
    [Displacement, checkDisp]=ReadPolyTech('Displacement.txt');
    [Current, checkCurr]=ReadPolyTech('Current.txt');
    Output=struct;
    if checkCurr==checkDisp
        disp(strcat('Matching Labels, ignoring measurement: ',freq,'-', volt));
    else
        
        Output.current=Current;
        Output.displacement=Displacement;
        Output.inputVoltage=volt;
        Output.AmplifierGain=2.1825066666;
    end
end


% Function designed to convert the txt format from polytechs output
% to matlab compatible data, in the form of a 2xn matrix. also outputs
% check to ensure that for every folder we are getting both current and
% displacement, not just two of either and 0 of the other.
function [Data, check]=ReadPolyTech(txtName)
% open file
    fid = fopen(txtName);
% pass through first line
    fgets(fid);
%     get second line
    check=fgets(fid);
%     get the first letter of measured channel, R for current, V for
%     displacement.
    check=check(16);
    %     pass through next 3 lines.
    fgets(fid);fgets(fid);fgets(fid);
%     retrieve all the numerical data
    Numbers=fscanf(fid,'%f');
%     close document
    fclose(fid);
%     convert 1x2n column two 2xn matrix, separating time and level data
    Data=vec2mat(Numbers,2);
    
end