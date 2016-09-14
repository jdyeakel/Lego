using Distributions
using Gadfly
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")

num_play = 20

init_probs = [
p_n=0.01,
p_a=0.01,
p_m=0.1,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]

int_m, sp_m, t_m, tp_m = build_template_degrees(num_play,init_probs);
writedlm("$(homedir())/Dropbox/PostDoc/2014_Lego/data_template/fw.csv",t_m);

keeprc=find(x->x>0,sum(t_m,1));


# Assessing stats across community size ranges
svecpower = collect(1:0.05:4.0);
svec=Array{Int64}(length(svecpower))*0;
for i=1:length(svecpower)
  svec[i]=round(Int64,10^svecpower[i]);
end
ls = length(svec);
C_out=zeros(ls);
init_probs = [
p_n=0.01,
p_a=0.01,
p_m=0.001,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
];
for i=1:ls
  s=svec[i];
  print(s);
  num_play = s;


  int_m, sp_m, t_m, tp_m = build_template_degrees(num_play,init_probs);
  #If there are no primary producers, then rinse and repeat
  # if sum(t_m[:,1])==0
  #   check = true;
  # else
  #   check = false;
  # end
  diag_int=diag(int_m);
  S = length(find(x->x=='n',diag_int));
  #Directed Connectance
  L = sum(t_m)/2;
  C = L/(S^2);
  C_out[i] = C;
end

sc_data = hcat(svec,C_out);
writedlm("$(homedir())/Dropbox/PostDoc/2014_Lego/data_template/size_conn.csv",sc_data);

sc_data=readdlm("$(homedir())/Dropbox/PostDoc/2014_Lego/data_template/size_conn.csv");

SCplot=plot(x=sc_data[:,1],y=sc_data[:,2],Geom.point,Scale.x_log10,Scale.y_log10,Theme(default_point_size=1pt, highlight_width = 0pt),Guide.xlabel("S"),Guide.ylabel("C"));

draw(PDF("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_SvC.pdf", 4inch, 3inch), SCplot)
