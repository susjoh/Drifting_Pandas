library(shiny)


##########################################################################
#
# Simulating a small population with drift and genotypic selection in a randomly mating dioecious species
# Written by Darren Obbard June 2020 to replace an EEG3 practical during COVID!
# Original basis of the Red Panda practical (probably) designed by Richard Ennos in the early 2000's (?)
#
##########################################################################

#this is the whole simulation for a single replicate, everything else is formatting and input/output
drift.selection.addrep<-function(p0=1,Ne=100,w=c(1,1,1),ngen=400){
    #always have at least one of the ancestral allele
    start_count<-min(p0,(2*Ne-1))
    #half the population will be each sex. Could alter the simulation to allow selfing or hermaphroditism
    m<-round(Ne/2,0)   
    #starting genotypes made by sampling alleles at random
    genotypes<-matrix(sample(c(rep(0,(2*Ne)-start_count),rep(1,start_count))),ncol=2)
    #vector to store mutant allele counts over time
    p<-rep(NA,ngen)
    for(i in 1:ngen){
        #store allele count for ith generation
        p[i]<-sum(genotypes)
        #calculate fitness
        fitness<-w[rowSums(genotypes)+1]
        #Red pandas cannot self!
        #mating is at random - half sibs are much more likely than full sibs
        #alleles are drawn with replacement from the parental population
        n_genotypes<-cbind(
            #sample maternal alleles (top half the popn are mothers)
            sample(c(genotypes[1:m,1],genotypes[1:m,2]),size=Ne,replace=TRUE, prob=c(fitness[1:m],fitness[1:m])),
            #sample paternal alleles  (bottom half the popn are father)
            sample(c(genotypes[(m+1):Ne,1],genotypes[(m+1):Ne,2]),size=Ne,replace=TRUE, prob=c(fitness[(m+1):Ne],fitness[(m+1):Ne]))
        )
        #replace the old genotypes with new genotypes
        genotypes<-n_genotypes
        #once we reach fixation or loss, no point in continuing the simulation. Save the current frequency in the last generation, and return the replicate
        if(p[i]==0|p[i]==(2*Ne)){
            p[ngen]<-p[i]
            return(p)
        }
    } #end gen
    #return the replicate
    return(p)
}


##########################################################################
#
# Shiny page layout
#
##########################################################################

ui <- fluidPage(
    
    titlePanel(
            div(img(height = 150,  src = "RedPanda.jpg"), 
               "Simulating a small population of Red Pandas")
            ),
    
    sidebarLayout(
        # panel with all inputs
        sidebarPanel(
            h4("Instructions:\n"),
            p("Follow the instructions in the practical guide to set up the simulation parameters, then press Run"),
            div(actionButton('go', 'Run'),align="right"),
         
            
            br(), 
            br(), 
            numericInput(inputId="Ne",label="Number of red pandas:",value=100,min=1,max=5000,width="210px"),
            numericInput(inputId="ngen",label="Number of generations:",value=200,min=1,max=10000,width="210px"),            
            numericInput(inputId="p0",label="Copies of the mutant allele:",value=100,min=1,width="210px"),
            numericInput(inputId="nrep",label="Number of simulations:",value=25,min=1,max=1000,width="210px"),
            p(strong("Relative genotype fitness (where the mutant allele is A):")),
            #fluidRow(
            splitLayout(
                numericInput(inputId="Waa",label="W[aa]",value=1.0,min=0,max=1,step=0.001,width="90px"),
                numericInput(inputId="WAa",label="W[Aa]",value=1.0,min=0,max=1,step=0.001,width="90px"),
                numericInput(inputId="WAA",label="W[AA]",value=0.98,min=0,max=1,step=0.001,width="90px")
            ),
            # ),
            br(), 
            br(), 
            
            # buttons to start, stop, reset
            p("Use these buttons to download the results and graphs for the simulation results shown"),
            fluidRow(
                column(5, downloadButton('graphs', 'Graphs')),
                column(5, downloadButton('results', 'Results'))
            ),                  
        ),
        
        # plot panel
        mainPanel(
            plotOutput("plot1",width="100%",height="650px"),
            hr(),
            fluidRow(
                textOutput("rep"),
                textOutput("lost"),
                textOutput("fixed"),
                textOutput("AlleleFreq")
            )

        )
    )
)


##########################################################################
#
# Shiny server functions
# Acheived through googling and copy&paste, so unlikely to be a good way of doing anything
#
##########################################################################

server <- function(input, output, session) {

    # reactive to store all reactive variables
    waits <- reactiveValues(rep=0, res=vector())
    
    # function to move simulation forward one replication
    forward <- function() {
            waits$rep <- waits$rep+1
            w<-isolate(c(input$Waa,input$WAa,input$WAA))
            w<-(w/max(w))
            p0<-isolate(input$p0)
            Ne<-isolate(input$Ne)
            ngen<-isolate(input$ngen)
            waits$res <-rbind(waits$res,drift.selection.addrep(p0=p0,Ne=Ne,w=w,ngen=ngen))
    }
    
    # when Run button is pressed
    observeEvent(input$go, {
		#slow the simulation down for small numbers of replicates, to allow the user to watch it.
        if(input$nrep<15){delay<-1000}
        if(input$nrep>=15&input$nrep<50){delay<-500}
        if(input$nrep>=50&input$nrep<=100){delay<-250}
        if(input$nrep>100){delay<-100}
        waits$rep<-0
        waits$res<-vector() 
        observe({
			#step the simulation forward
            isolate({
                forward()
            })
			# all the time there are still steps to run, ensure that the figures will be re-rendered after a delay
            if (isolate(waits$rep) < isolate(input$nrep)){
                invalidateLater(delay, session)
            }
        })
        waits$rep<-0
    })

    # render the plots, and generate the output numbers        
     output$plot1<-renderPlot({
        # the isolates()'s seem to prevent replotting on input changes, allowing the 'Run' button complete control
        Ne<-isolate(input$Ne)
        ngen<-isolate(input$ngen)
        nrep<-isolate(input$nrep)
        res<- waits$res
        #text output of replicate number
        output$rep<-renderText({paste("Replicate number ",waits$rep,sep="")})
        
        #render empty plots prior to running anything
        par(mar=c(6,5,4,4),mfrow=c(2,1),cex.axis=1.0, cex.lab=1.2, cex.main=1.5) 
        plot(0,0,type="n",xlim=c(1,ngen),ylim=c(0,1),xlab="Generations",ylab="Frequency of allele A",main="Allele frequency over time")
        barplot(rep(0,min(20,2*Ne)+2),col="gray30",ylab="Number of replicates",las=2,main="Histogram of final frequencies for allele A",space=0,ylim=c(0,1))
        all_cols<-rainbow(nrep)
        
        #render plots
        if(is.matrix(res)){
            plot(0,0,type="n",xlim=c(1,ngen),ylim=c(0,1),xlab="Generations",ylab="Frequency of allele A",main="Allele frequency over time")
            # if more than 100 replicates only plot the last 100 to reduce plotting time
            if(nrow(res)<=100){
                for(i in 1:nrow(res)){points(1:length(res[i,]),res[i,]/(2*Ne),type="l",lwd=2,col=adjustcolor(all_cols[i], alpha.f = 0.4))}
            }else{
                #plot the last 50, plus up to 50 non-losses before that.
                
                nl<-0
                for(i in (nrow(res)-50):1){
                                              if(res[i,ncol(res)]!=0 & nl<=50){
                                                 points(1:length(res[i,]),res[i,]/(2*Ne),type="l",lwd=2,col=adjustcolor(all_cols[i], alpha.f = 0.4))
                                                  nl<-nl+1
                                              }
                }
                for(i in (nrow(res)-50):nrow(res)){points(1:length(res[i,]),res[i,]/(2*Ne),type="l",lwd=2,col=adjustcolor(all_cols[i], alpha.f = 0.4))}

                
            }
            
            #Different barchart intervals for small population and large population sizes
            if(Ne<10){
                #small sample, use direct counts
                counts<-table(c(res[,ncol(res)],0:(2*Ne)))-1
                names(counts)<-round(as.numeric(names(counts))/(2*Ne),2);names(counts)[1]<-"lost";names(counts)[length(counts)]<-"fixed"
                mids<-barplot(counts,col=c("red",rep("gray30",length(counts)-2),"red"),ylab="Number of replicates",las=2,main="Histogram of final frequencies for allele A",space=0)
                #text output coutns of loses and fixations
                output$lost<-renderText({paste("Allele A has been LOST ",100*signif(counts[1]/nrow(res),3),"% of the time",sep="")}) 
                output$fixed<-renderText({paste("Allele A has been FIXED ",100*signif(counts[length(counts)]/nrow(res),3),"% of the time",sep="")})               
                
            }else{
                #large sample, use hist (without plotting it) to carve up the counts into buckets
                # with extra effort to make the first and last buckets count losses and fixations as special cases
                breaks<-seq(from=0,to=(2*Ne),length=20);breaks[1]<-0.5;breaks[length(breaks)]<-(2*Ne)-0.5;
                breaks<-c(0,breaks,2*Ne)
                y<-hist(res[,ncol(res)],breaks=breaks,plot=FALSE)$counts
                labels<-paste("<",as.character(round(breaks[-1]/(2*Ne),2)),sep="");labels[1]<-"lost";labels[length(labels)]<-"fixed"
                mids<-barplot(y,col=c("red",rep("gray30",length(y)-2),"red"),ylab="Number of replicates",las=2,names.arg=labels,main="Histogram of final frequencies for allele A",space=0)
                #text output coutns of loses and fixations
                output$lost<-renderText({paste("Allele A has been LOST ",100*signif(y[1]/nrow(res),3),"% of the time, i.e. on ",y[1]," occasions",sep="") })
                output$fixed<-renderText({paste("Allele A has been FIXED ",100*signif(y[length(y)]/nrow(res),3),"% of the time, i.e. on ",y[length(y)]," occasions",sep="")  })                
            }
          
            freqs<-res[,ncol(res)]/(2*Ne)
            # Calculating diversity is now an excercise for the reader
            # diversites<-1-  (freqs^2 + (1-freqs)^2)
            
            # text output for the screen
            output$AlleleFreq<-renderText({paste("Mean frequency of A across replicate end-points was ",round(100*mean(freqs),2),"%",sep="")})    

        }
       
    }) 
     
     #download the frequencies from the simulation
     output$results <- downloadHandler(
         filename = function() { paste("AlleleFrequenciesOverTime.",input$Ne,"_Pandas.",input$p0,"_StartingMutations.GenotypeFitnesses_",paste(c(input$Waa,input$WAa,input$WAA),collapse="-"),'.tsv', sep='') },
         content = function(file) {
             res<-waits$res
             if(is.matrix(res)){
                rownames(res)<-paste("rep",1:nrow(res),sep="_")
                colnames(res)<-paste("gen",1:ncol(res),sep="_")
             }
             write.table(res/(2*input$Ne), file,sep="\t",quote=FALSE)
         }
     )
     
     #download the plots, after re-ploting as a png
     output$graphs <- downloadHandler(
         filename = function() { paste("Graphs.",input$Ne,"_Pandas.",input$p0,"_StartingMutations.GenotypeFitnesses_",paste(c(input$Waa,input$WAa,input$WAA),collapse="-"),'.png', sep='') },
         content = function(file) {
             #set up the png size and format  
             png(file,  width= 4,height = 6,units="in", res= 300, pointsize = 7)
             #replot the graphs exactly as before
             Ne<-isolate(input$Ne)
             ngen<-isolate(input$ngen)
             nrep<-isolate(input$nrep)
             res<- waits$res     
             par(mar=c(6,5,4,4),mfrow=c(2,1),cex.axis=1.0, cex.lab=1.2, cex.main=1.5) 
             all_cols<-rainbow(nrep)
             if(is.matrix(res)){
                 plot(0,0,type="n",xlim=c(1,ngen),ylim=c(0,1),xlab="Generations",ylab="Frequency of allele A",main="Allele frequency over time")
                 if(nrow(res)<=100){
                     for(i in 1:nrow(res)){points(1:length(res[i,]),res[i,]/(2*Ne),type="l",lwd=2,col=adjustcolor(all_cols[i], alpha.f = 0.4))}
                 }else{
                     for(i in (nrow(res)-100):nrow(res)){points(1:length(res[i,]),res[i,]/(2*Ne),type="l",lwd=2,col=adjustcolor(all_cols[i], alpha.f = 0.4))}
                 }
                 if(Ne<10){ 
                     counts<-table(c(res[,ncol(res)],0:(2*Ne)))-1
                     names(counts)<-round(as.numeric(names(counts))/(2*Ne),2);names(counts)[1]<-"lost";names(counts)[length(counts)]<-"fixed"
                     mids<-barplot(counts,col=c("red",rep("gray30",length(counts)-2),"red"),ylab="Number of replicates",las=2,main="Histogram of final frequencies for allele A",space=0)
                 }else{
                     breaks<-seq(from=0,to=(2*Ne),length=20);breaks[1]<-0.5;breaks[length(breaks)]<-(2*Ne)-0.5;
                     breaks<-c(0,breaks,2*Ne)
                     y<-hist(res[,ncol(res)],breaks=breaks,plot=FALSE)$counts
                     labels<-paste("<",as.character(round(breaks[-1]/(2*Ne),2)),sep="");labels[1]<-"lost";labels[length(labels)]<-"fixed"
                     mids<-barplot(y,col=c("red",rep("gray30",length(y)-2),"red"),ylab="Number of replicates",las=2,names.arg=labels,main="Histogram of final frequencies for allele A",space=0)
                 }
             }             
            #turn of the png             
            dev.off()
         }
    ) 
    
  
}

#run the server
shinyApp(ui, server)

