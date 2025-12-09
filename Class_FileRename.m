classdef Class_FileRename
    %CLASS_FIRERENAME Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        TopFolder  % Data folder with date
        ImageFolder
        PressureFolder
        PhotodiodeFolder

        ImageFolderPattern = 'Image';
        PressureFolderPattern = 'Pressure';
        PhotodiodeFolderPattern = 'Photodiode';
        
        bRenamePressure = 1;
        bRenamePhotodiode = 1;

        ImageFiles

        PressureFiles
        PressureIndexFiles
        PressureLogFiles

        PhotodiodeFiles
        PhotodiodeIndexFiles
    end
    
    methods
        function obj = Class_FileRename(topFolder)
            if ~isfolder(topFolder)
                error('Input folder name is not a folder')
            end
            obj.TopFolder = topFolder;
            obj = obj.FindSubFolder();
            obj = obj.CompareFileAmount();
            obj = obj.RenameFiles();
        end
        %%
        function obj = FindSubFolder(obj)
            topFolder = obj.TopFolder;
            SubDir = dir(topFolder);
            SubDirName = {SubDir.name};
            for ii = 1:numel(SubDirName)
                if contains(SubDirName{ii}, obj.PressureFolderPattern)
                    obj.PressureFolder = fullfile(topFolder, SubDirName{ii});
                end
                if contains(SubDirName{ii}, obj.PhotodiodeFolderPattern)
                    obj.PhotodiodeFolder = fullfile(topFolder, SubDirName{ii});
                end
                if contains(SubDirName{ii}, obj.ImageFolderPattern)
                    obj.ImageFolder = fullfile(topFolder, SubDirName{ii});
                end
            end

            obj = obj.GetImageFileList();
            obj = obj.GetPressureFileList();
            obj = obj.GetPhotodiodeFileList();
        end
        %%
        function obj = GetImageFileList(obj)
            if isempty(obj.ImageFolder)
                warning('Image folder not found in: %s', obj.TopFolder)
                return
            end            
            obj.ImageFiles = dir(fullfile(obj.ImageFolder,'*.cine'));
            
        end

        function obj = GetPressureFileList(obj)
            if isempty(obj.PressureFolder)
                warning('Pressure folder not found in: %s', obj.TopFolder)
                return
            end            
            obj.PressureFiles = dir(fullfile(obj.PressureFolder,'*.tdms'));
            obj.PressureIndexFiles = dir(fullfile(obj.PressureFolder,'*.tdms_index'));
            obj.PressureLogFiles = dir(fullfile(obj.PressureFolder,'*.log'));
        end

        function obj = GetPhotodiodeFileList(obj)
            if isempty(obj.PhotodiodeFolder)
                warning('Photodiode folder not found in: %s', obj.TopFolder)
                return
            end            
            obj.PhotodiodeFiles = dir(fullfile(obj.PhotodiodeFolder,'*.tdms'));
            obj.PhotodiodeIndexFiles = dir(fullfile(obj.PhotodiodeFolder,'*.tdms_index'));
        end
    
        function obj = CompareFileAmount(obj)
            if length(obj.PressureFiles) ~= length(obj.PressureIndexFiles) && obj.bRenamePressure == 1
                warning('.tdms file amount != .tdms_index file amount, abort renaming')
                return
            end

            if length(obj.PressureFiles) ~= length(obj.ImageFiles) && obj.bRenamePressure == 1
                warning('.tdms file amount != .cine file amount, abort renaming')
                return
            end

            if length(obj.PhotodiodeFiles) ~= length(obj.ImageFiles) && obj.bRenamePhotodiode == 1
                warning('.tdms (Photodiode) file amount != .cine file amount, abort renaming')
                return
            end

            fprintf('File names match, can proceed with renaming\n')
        end
%%
        function obj = RenameFiles(obj)
            ImageFileNames = {obj.ImageFiles.name}; 
            for jj = 1:length(ImageFileNames)
                ImageFileNames{jj} = strrep(ImageFileNames{jj}, '.cine','');
            end

            filesToRename_Pressure = {obj.PressureFiles.name}; 
            filesToRename_Pressure_Index = {obj.PressureIndexFiles.name}; 
            filesToRename_Pressure_log  = {obj.PressureLogFiles.name}; 
            
            if obj.bRenamePressure == 1
                for ii = 1:length(filesToRename_Pressure)
                    % Extract the directory and name from the files to rename
                    [pathstr,~,~] = fileparts(fullfile(obj.PressureFiles(1).folder,filesToRename_Pressure{ii}));
                    
                    % Construct the new file name with the desired extension
                    newFileName = fullfile(pathstr, [ImageFileNames{ii}, '.tdms']);
                    if isfile(newFileName)
                        fprintf('File exist: %s\n', newFileName);
                        continue
                    end
                    % Rename the file
                    movefile(fullfile(obj.PressureFiles(1).folder,filesToRename_Pressure{ii}), newFileName);
                end
    
                for ii = 1:length(filesToRename_Pressure_Index) 
                    % Extract the directory and name from the files to rename
                    [pathstr,~,~] = fileparts(fullfile(obj.PressureFiles(1).folder,filesToRename_Pressure_Index{ii}));
                    
                    % Construct the new file name with the desired extension
                    newFileName = fullfile(pathstr, [ImageFileNames{ii}, '.tdms_index']);
                    if isfile(newFileName)
                        fprintf('File exist: %s\n', newFileName);
                        continue
                    end
                    % Rename the file
                    movefile(fullfile(obj.PressureFiles(1).folder,filesToRename_Pressure_Index{ii}), newFileName);
                end
    
                for ii = 1:length(filesToRename_Pressure_log)  
                    % Extract the directory and name from the files to rename
                    [pathstr,~,~] = fileparts(fullfile(obj.PressureFiles(1).folder,filesToRename_Pressure_log{ii}));
                    
                    % Construct the new file name with the desired extension
                    newFileName = fullfile(pathstr, [ImageFileNames{ii}, '.log']);
                    if isfile(newFileName)
                        fprintf('File exist: %s\n', newFileName);
                        continue
                    end
                    % Rename the file
                    movefile(fullfile(obj.PressureFiles(1).folder,filesToRename_Pressure_log{ii}), newFileName);
                end         
            end


            % Photodiode
            filesToRename_Photodiode  = {obj.PhotodiodeFiles.name}; 
            filesToRename_Photodiode_Index = {obj.PhotodiodeIndexFiles.name}; 
            if obj.bRenamePhotodiode == 1
                for ii = 1:length(filesToRename_Photodiode)  
                    % Extract the directory and name from the files to rename
                    [pathstr,~,~] = fileparts(fullfile(obj.PhotodiodeFiles(1).folder,filesToRename_Photodiode{ii}));
                    
                    % Construct the new file name with the desired extension
                    newFileName = fullfile(pathstr, [ImageFileNames{ii}, '.tdms']);
                    if isfile(newFileName)
                        fprintf('File exist: %s\n', newFileName);
                        continue
                    end
                    % Rename the file
                    movefile(fullfile(obj.PhotodiodeFiles(1).folder,filesToRename_Photodiode{ii}), newFileName);
                end   
            
                for ii = 1:length(filesToRename_Photodiode) 
                    % Extract the directory and name from the files to rename
                    [pathstr,~,~] = fileparts(fullfile(obj.PhotodiodeFiles(1).folder,filesToRename_Photodiode{ii}));
                    
                    % Construct the new file name with the desired extension
                    newFileName = fullfile(pathstr, [ImageFileNames{ii}, '.tdms_index']);
                    if isfile(newFileName)
                        fprintf('File exist: %s\n', newFileName);
                        continue
                    end
                    % Rename the file
                    movefile(fullfile(obj.PhotodiodeFiles(1).folder,filesToRename_Photodiode_Index{ii}), newFileName);
                end
            end
        end
    end
end



