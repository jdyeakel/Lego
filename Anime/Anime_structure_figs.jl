loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

d = load(string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/structure.jld"));
#This loads the dictionary
rich = d["rich"];
sprich = d["sprich"];
turnover = d["turnover"];
res_overlap = d["res_overlap"];
conn = d["conn"];
prim_ext = d["prim_ext"];
sec_ext = d["sec_ext"];
status = d["status"];
lpot_col = d["lpot_col"];
degrees = d["degrees"];
trophic = d["trophic"];


