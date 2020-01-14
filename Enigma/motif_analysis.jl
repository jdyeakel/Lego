using RCall

M4_1 = [0 0 0 0;0 0 0 0;0 0 1 1;0 0 1 0];

M4_2 = [0 0 0 0;0 0 0 0;0 0 1 0;0 0 1 0];

T4_1 = [0 0 0 0;0 0 0 0;0 0 1 1;0 0 1 0];

T4_2 = [0 0 0 0;0 0 0 0;0 0 0 0;0 0 1 0];


M6_1 = [0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 1 1 1;0 0 0 1 1 0;0 0 0 1 0 0];

M6_2 = [0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 1 1 1;0 0 0 1 1 0;0 0 0 0 0 0];

T6_1 = [0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 1 1 1;0 0 0 1 1 0;0 0 0 1 0 0];

T6_2 = [0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 1 1 0;0 0 0 1 0 0];

motiflist = [M6_1,M6_2,T6_1,T6_2];

motiflist = [M4_1,M4_2,T4_1,T4_2];

nestedness = Array{Float64}(undef,4);
for i=1:4
    R"""
    library(UNODF)
    library(igraph)
    unodfvalue = suppressWarnings(unodf($(motiflist[i]),selfloop=TRUE))
    comb_nestvalue_c = unodfvalue$UNODFc
    comb_nestvalue_r = unodfvalue$UNODFr
    """
    @rget comb_nestvalue_c;
    @rget comb_nestvalue_r;
    nestedness[i] = mean([comb_nestvalue_c,comb_nestvalue_r])
end
    
R"""
plot(c(1,2),$(nestedness[1:2]),xlim=c(0.5,2.5),ylim=c(0,0.5),xlab='Time step',ylab='Nestedness (UNODF)')
lines(c(1,2),$(nestedness[1:2]))
points(c(1,2),$(nestedness[3:4]),pch=2)
lines(c(1,2),$(nestedness[3:4]),lty=2)
"""

R"""
image($M6_1[1:6,6:1],xaxt='n',yaxt='n',col=c('white','black'))
"""

R"""
image($M6_2[1:6,6:1],xaxt='n',yaxt='n',col=c('white','black'))
"""

R"""
image($T6_1[1:6,6:1],xaxt='n',yaxt='n',col=c('white','black'))
"""

R"""
image($T6_2[1:6,6:1],xaxt='n',yaxt='n',col=c('white','black'))
"""

