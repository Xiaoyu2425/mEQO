# Ensemble Quotient Optimization for microbiome <img src="https://tva1.sinaimg.cn/large/008vxvgGgy1h7di5wet89j30e00g8mxp.jpg" align="right" height="137" />

**"Harmony without uniformity."**   --*The Analects* ~ 400 B.C



## 1. Introduction

Unsupervised microrbiome coarse-graining can be generalized using the optimization of Ensemble Quotient. This R package *mEQ*O enables microbiome Ensemble Quotient Optimization (EQO) guided by trait variables, which can be continuous, discrete, or even uniform. 

In short, EQO finds an assemblage of species which as a whole is most strongly indicative to a trait variable, even if the individual species may be decoupled from the trait variable. This marks a fundamental distinction between EQO and previous clustering algorithms. Previous clustering algorithms usually search for a group of species that are intercorrelated with each other, while EQO finds species that are **"consistent but complement each other"**, so that they as a whole are strongly coupled with a trait variable. Full details about the mathematical framework, statistical interpretations and algorithmic design can be found in this [manuscript](https://www.biorxiv.org/content/10.1101/2022.08.02.502537v1).

## 2. Installation

*mEQO* can be installed by 

```R
devtools::install_github("Xiaoyu2425/mEQO")
```

Note: These R packages will be automatically imported as dependencies of mEQO: GA, doParallel, data.table, reshape2. In addition, if you want to implement EQO with the boolean least square (BLS) algorithm, a Gurobi optimizer and its R interface will be required as dependencies. A free license for academic use can be easily acquired on their [website](https://www.gurobi.com). An installation guide can be found [here](https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation_guide.html). 

## 3. Example 

Input for EQO is simply the microbiome compositon and the functional/environmental variable. Microbiome composition should be a matrix whose rows are samples and columns are taxa. Column names should be provided. Functional/environmental variable should be a vector whose length is the number of samples. You can find an example (`Microbiome` and `trait`) for input format after loading mEQO by running `library(mEQO)`.

Here we demonstrate the use of mEQO with this example of simplified dataset. As shown in the figure below, in this mock microbiome, the relative abundance of single species (color strips) is uncorrelated with a functional trait (black dots). The goal of the algorithm is to identify a group of species whose total abundance is most strongly correlated with the functional trait. 

<img src="https://tva1.sinaimg.cn/large/008vxvgGgy1h7dniv5szcj31900u0whr.jpg" alt="composition1" height="300" width="450" />

To that end, we can run one of the functions below to implement EQO.  

```R
EQO_bls(Microbiome,trait) # implemented by BLS
EQO_ga("c",Microbiome,trait) # or, implemented by GA
```

We have provided two different algorithms for implementing EQO, each with its advantages and disadvantages: boolean least square (BLS) and genetic algorithm (GA). BLS works very accurately for smaller-scale problems with a continuous variable (e.g., < 300 species). GA can be more efficient for larger-scale datasets, whose accuracy depends on number of iterations. Please refer to "additional discussion" for full details. 

In the above example of GA, "c" indicates a continuous trait variable. In addition, `EQO_ga()` can also support a trait variable that is discrete or uniform. You can specify "d" for a discrete trait variable (e.g., healthy or diseased samples) or "u" for a uniform trait variable (e.g., stable community composition across samples).

Here, the best group of species found by the algorithm that is coupled to the trait variable is comprised of species 1, 2 and 3, as shown in the figure below. 

<img src="https://tva1.sinaimg.cn/large/008vxvgGgy1h7dniyrywvj31900u0dhz.jpg" alt="composition2" height="300" width="450" />

For microbiome studies, in addition to obtaining a group of species with strongest statistical power, we may be more interested in understanding which single species or what combinations of species within that group have the highest importance. 

We have two suggestions in this regard. 

Firstly, you can apply regularization to the algorithm by specifying *Nmax* as the maximal number of species allowed in the group. In statistics, this is often based on Akaike Information Criterion (AIC). Here, the minimal AIC is achieved at *Nmax=3*, as shown in the figure below. 

```R
aic<-sapply(1:5,function(N){
	assemblage<-EQO_bls("c",Microbiome,trait,Nmax=N)$abundance
	return(AIC(lm(trait~assemblage))+2*(N-1))
})
```

<img src="https://tva1.sinaimg.cn/large/008vxvgGgy1h7dexjolgcj31900u0wfy.jpg" alt="aic" height="280" width="400" align="center" />

Secondly, we recommend illustrating your functional group with a Cross-validation-based Aggregation Network (CAN). We provide a `CAN()` function in the *mEQO* package for that. In short, this function performs the following: (1) splitting your original dataset into a test set and a validation set, (2) performing EQO with the test set and cross-validation with the validation set, 3) evaluating relative importance of individual species or pairs by repeating (1) and (2) for multiple times. Remarkably, CAN can illustrate patterns that are very different from a traditional association network such as a co-occurrence network. In CAN, each nodes still denotes a species, while the weight of an edge indicates the strength of cross-validated correlation with the functional trait when the two connected species arse coarse-grained in a group. As is shown in the following figure, species 1, 2, and 3 stand out from the cross-validation as a highly interconnected module, strongly suggesting their emergent ecological role as a group together. 

```R
can<-CAN("c",Microbiome,trait)

nodes<-can$nodes
nodes$color.background=c("#CF0A0A","#F8C957","#125D98","#928B8B","#C65D7B","#4E944F","#5EA3A6","#DD6B4D")
nodes$color.border <- "black"

edges<-can$edges
edges$color<-"#99C4C8"

library(visNetwork)
visNetwork(nodes,edges) %>%
	visConfigure(enabled = TRUE) %>% 
	visLayout(randomSeed = 123) %>%
	visPhysics(enabled=F)
```



<img src="https://tva1.sinaimg.cn/large/008vxvgGgy1h7dnnvqm9hj31520lqq5d.jpg" alt="network2" height="280" width="500" />



## 4. Additional discussion

The major computational challenge for microbiome coarse-graining lies in its combinatorial optimization nature, which is still one of the major open challenges in the modern operational research. In the *mEQO* package, we have tackled this challgenge with BLS and GA as two major approaches. 

For the BLS algorithm, we leverage a mature commercial mixed interger programming solver, the Gurobi optimizer, to implement EQO. Please note that BLS only works for a continuous trait variable and relatively small-scale problems (e.g., < 300 species). In order to run BLS more efficiently, we strongly recommend reducing the total number of species, e.g., by removing the long tail of singleton species. Applying regularization to the maximal number of species allowed is also very important to make BLS more efficient. 

GA can be more suitable for larger-scale problems. [GA](https://cran.r-project.org/web/packages/GA/vignettes/GA.html) has been widely used as a successful strategy for stochastic binary searching with convergence in probability. However, the rate of convergence is dependent on size and structure of spefific datasets. Below is an illustration of the accuracy of GA at different number of iterations with simulated microbiome datasets. 

<img src="https://tva1.sinaimg.cn/large/008vxvgGgy1h7hawyk0hrj30pq0roq4r.jpg" alt="??????2022-10-24 ??????9.59.06" height="400" width="380" />

You can determine the optimal number of iterations for your specific datasets based on your computational power and desirable accuracy, for example

```R
iterations<-c(100,200,500,1000)
fitness<-replicate(10,sapply(iterations,function(t){EQO_ga("c",Microbiome,trait,maxIter=t)$fitness}))
plot(iterations,apply(fitness,2,mean))
```

Importantly, Instead of focusing on absolute optimality (e.g., having computing time doubled to optimize from r = 0.80 to r = 0.81),  it is of greater ecological significance to ask what species combinations are the most robust and important across multiple instances of cross-validations. That is also why we recommend constructing a CAN with our approach, especially for a large-scale dataset for which combinatorial optimization struggles to rapidly reach an exact optimal solution.



## 5. Partially known functional groups

mEQO also supports cases where you want to designate certain species to be included in the functional group, when you have *a priori* knowledge of the functional role of that species. Then mEQO can start on the basis of that partially known functional group and continue optimizing by combining other species. In the package, this can be implemented by specifying the argument *pk* for function `EQO_ga()`. As a vector whose length is equal to the total number of species in the microbiome, *pk* indicates partially known functional group based on a priori knowledge, with 1 for species forced to be included in the targeted group and 0 for the other unknown species. EQO will then search on the basis of the provided partially known group without removing the designated species. For instance, in the case where we want to ensure species 4 is included in the final group, we would include the following argument:

```R
EQO_ga("c",Microbiome,trait,pk=c(0,0,0,1,0,0,0,0))
```



## 6. Citation

Please cite [mEQO](https://www.biorxiv.org/content/10.1101/2022.08.02.502537v1.abstract) as an R package, as well as the following dependencies: [GA](https://cran.r-project.org/web/packages/GA/index.html); [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html); [data.table](https://cran.r-project.org/web/packages/data.table/index.html); [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html).

