function data = load_raw(file_name)

bytes_samp = 2; %How many bytes per sample for the int16  type

file_info = dir(file_name); %Info about the file
nSampsTotal = file_info.bytes/bytes_samp; %Find the total number of samples in the file to be analyzed

fid = fopen(file_name, 'r'); %Open an fid for the file of interest. This is something like a pointer for this file
data = fread(fid, [nSampsTotal,1],'int16'); %Get the data for the specified time period from the bin file for all channels

fclose(fid);



