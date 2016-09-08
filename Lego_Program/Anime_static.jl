using Distributions
include("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")

num_play = 100

init_probs = [
p_n=0.05,
p_a=0.02,
p_m=0.001,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]

int_m=build_template_degrees(num_play,init_probs);

#Derive the trophic adjacency matrix, which includes all 'ai','aa','an' interactions
t_m = zeros(num_play,num_play);
for i=1:num_play
  for j=1:num_play
    if i > j
      if int_m[i,j] == 'a' && int_m[j,i] == 'i'
        t_m[i,j] = 1
        t_m[j,i] = 1
      end
      if int_m[i,j] == 'i' && int_m[j,i] == 'a'
        t_m[i,j] = 1
        t_m[j,i] = 1
      end
      if int_m[i,j] == 'a' && int_m[j,i] == 'n'
        t_m[i,j] = 1
        t_m[j,i] = 1
      end
      if int_m[i,j] == 'n' && int_m[j,i] == 'a'
        t_m[i,j] = 1
        t_m[i,j] = 1
      end
      if int_m[i,j] == 'a' && int_m[j,i] == 'a'
        t_m[i,j] = 1
        t_m[j,i] = 1
      end
    end
  end
end

#Food web network
trophic_m = Array(Int,num_play,num_play)*0;
trophic = find(x->x=='a',int_m);
trophic_m[trophic]=1;
trophic_m
