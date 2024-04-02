function new_struct = concatStructs( struct1, struct2 )
%CONCATSTRUCTS makes one structure with the fields that existed in either
%structure. If the same field exists in both structures, it will be written
%over by the field in struct2.
new_struct = struct1;
struct2_fields = fields(struct2);
for i=1:length(struct2_fields)
    new_struct.(struct2_fields{i}) = struct2.(struct2_fields{i});
end

