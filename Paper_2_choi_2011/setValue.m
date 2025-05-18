function [] = setValue(field, var, val)

filename = "dimData.xml";

data = readstruct(filename);

data.(field).(var) = val;

writestruct(data, filename);

end