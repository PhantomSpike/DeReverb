function adj_coeff = db_adjust(data,db_new)

p0 = 20.*10.^-6; %p0 is the reference pressure in Pa used to calculate the dB
rms_old = rms(data);
rms_new = p0.*10.^(db_new/20);
adj_coeff = rms_new./rms_old;

end