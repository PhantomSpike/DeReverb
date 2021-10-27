function db = db_calc(data)

rms_data = rms(data);
p0 = 20.*10.^-6; %p0 is the reference pressure in Pa used to calculate the dB
db = 20.*log10(rms_data./p0);

end