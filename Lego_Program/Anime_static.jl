using Distributions
using Gadfly
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")

num_play = 100

init_probs = [
p_n=0,
p_a=0.01,
p_m=0,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]

int_m, t_m = build_template_degrees(num_play,init_probs);
writedlm("$(homedir())/Dropbox/PostDoc/2014_Lego/data_template/fw.csv",t_m);



# Assessing stats across community size ranges
svecpower = collect(1:0.02:4);
svec=Array{Int64}(length(svecpower))*0;
for i=1:length(svecpower)
  svec[i]=round(Int64,10^svecpower[i]);
end
ls = length(svec);
C_out=zeros(ls);
for i=1:ls
  s=svec[i];
  print(s);
  check=true;
  while check==true
    num_play = s;
    init_probs = [
    p_n=0.0,
    p_a=0.01,
    p_m=0.0,
    p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
    ];

    int_m, t_m = build_template_degrees(num_play,init_probs);
    #If there are no primary producers, then rinse and repeat
    # if sum(t_m[:,1])==0
    #   check = true;
    # else
    #   check = false;
    # end
  end
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

SCplot=plot(x=sc_data[:,1],y=sc_data[:,2],Geom.point,Scale.x_log10,Theme(default_point_size=1pt, highlight_width = 0pt),Guide.xlabel("S"),Guide.ylabel("C"));

draw(PDF("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_SvC.pdf", 4inch, 3inch), SCplot)
