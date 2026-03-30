plot_domrev <-function(h1,s1,h2,s2,maxfreq,reftime){
    
    # maxfreq -- max frequecy for y axis plotting
    # reference generation to measure prob established (i.e. prob established after x gen)-- cannnot exceed ngen
    
    N<-1000000                    # population size
    n_fix<-0                  # track fixations
    n_loss<-0                 # track losses
    q_eqb<-0                   # eqb frequency
    ngen<-1000                # number of generations
    nreps<-10000               # replicates of evolution
    
    # Relative fitnesses based on hs here
    #s1=(0)           # enter selection coefficient             DOMINANT ADVANTAGE
    #h1=1               # enter dominance
    w_AA = 1            # fitness of ancestral homozygote
    w_Aa = (1+h1*s1)      # fitness of heterozygote
    w_aa = (1+s1)        # fitness of mutant homozygote
    
    #s2=(-1)           # enter selection coefficient             RECESSIVE DELETERIOUS
    #h2=0               # enter dominance
    w_BB = 1            # fitness of ancestral homozygote
    w_Bb = (1+h2*s2)      # fitness of heterozygote
    w_bb = (1+s2)        # fitness of mutant homozygote
    
    
    for (j in 1:nreps){                # replicates of evolution
        q = 1/(2*N)                     # initial frequency
        freqs<-q                    # initialize a vector of frequencies
        for (i in 1:ngen){          # run each replicate for ngen
            
            # selection part
            mean_w = (((1-q)^2)*(w_AA)) + 2*q*(1-q)*(w_Aa) + q^2*(w_aa)     # calculate mean fitness
            q = (q*(1-q)*(w_Aa) + q^2*( w_aa))/mean_w                       # calculate new frequency
            
            mean_w = (((1-q)^2)*(w_BB)) + 2*q*(1-q)*(w_Bb) + q^2*(w_bb)     # calculate mean fitness
            q = (q*(1-q)*(w_Bb) + q^2*( w_bb))/mean_w                       # calculate new frequency
            
            # drift part
            A<-(rbinom(1, 2*N, q))  # number of A alleles
            q<-A/(2*N)              # frequency of A
            freqs<-c(freqs, q)
           
           if (i == reftime){ # log number of losses within reference time
               if (q==1){n_fix=n_fix+1}
               if (q==0){n_loss=n_loss+1}
               }
        } # ngen loop
        
        
        
        # plot
        if (j==1){plot(as.vector(freqs), type="l", ylim=c(0,maxfreq), ylab="frequency of a", xlab="generation", main=paste("Dominance-Reversal selection with drift, N = ", N, "\ndominant advantage hs = ", h1*s1, " recessive detriment h= ", h2, ", s=", s2))}
        else{lines(as.vector(freqs), type="l", ylim=c(0,1), col=j)}
        
        if (q !=0){q_eqb=q_eqb+q};
        
    } # replicates of evolution loop
    text(200,maxfreq-0.2*maxfreq, paste("prob_estab=", (1-(n_loss/nreps)), "\neqb freq=", round(q_eqb/(nreps-n_loss),3)))
    
}


