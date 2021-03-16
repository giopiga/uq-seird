%% In order to correctly install all the folders, you need to add all the paths
% run this script to set them (already run in the example_main.m)
clear
close all
clc

my_dir = fileparts(which("SEIRDcode_start.m")); 
addpath(genpath(my_dir));
rmpath(strcat(my_dir,"\Code\InvUQ\MCMCstat"));

FwdUQ_dir=strcat(my_dir,"\Code\FwdUQ");
InvUQ_dir=strcat(my_dir,"\Code\InvUQ");
Models_dir=strcat(my_dir,"\Code\Models");
SA_dir=strcat(my_dir,"\Code\SA");
Utils_dir=strcat(my_dir,"\Code\Utils");

addedPath=0;
addedPath=addedPath+(~isempty(strfind(path,FwdUQ_dir)));
addedPath=addedPath+(~isempty(strfind(path,InvUQ_dir)));
addedPath=addedPath+(~isempty(strfind(path,Models_dir)));
addedPath=addedPath+(~isempty(strfind(path,SA_dir)));
addedPath=addedPath+(~isempty(strfind(path,Utils_dir)));

if addedPath==5
    disp("Everything has been set correctly")
    pause(0.5)
    clc
    disp("Everything has been set correctly.")
    pause(0.5)
    clc
    disp("Everything has been set correctly..")
    pause(0.5)
    clc
    disp("Everything has been set correctly...")
    pause(0.5)
    
    clc    
else
    error("Ouch... Something went wrong with the paths configuration: try again or set them by hands.")
end

clear
close all