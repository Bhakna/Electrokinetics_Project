function [] = setValue(datatype, field, var, val)

if datatype == "dim"
    filename = "dimData.xml";
    data = readstruct(filename);
    data.(field).(var) = val;
    writestruct(data, filename);
    nondimData();
elseif datatype == "nondim"
    filename = "nondimData.xml";
    data = readstruct(filename);
    data.(var) = val;
    writestruct(data, filename);
elseif datatype == "default"
    filename = "dimData.xml";
    data = readstruct("defaultData.xml");
    writestruct(data, filename);
    nondimData();
else
    error("Incorrect DataType. Please use 'dim', 'nondim', or 'default' ");
end

end