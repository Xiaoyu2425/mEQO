# Ensemble Quotient Optimization for microbiome (mEQO) <img src="https://tva1.sinaimg.cn/large/008vxvgGgy1h7di5wet89j30e00g8mxp.jpg" align="right" height="137" />

**"Harmony without uniformity."** --*The Analects* ~ 400 B.C

## 1. Introduction

EQO finds an assemblage of species which as a whole is most strongly indicative to a trait variable, even if the individual species therein is not. In contrast to most clustering algorithms that usually search for a group of species intercorrelated with each other, EQO finds species that are **"consistent but complement each other"**, so that they as a whole are strongly coupled with a trait variable. Full details about the mathematical framework, statistical interpretations and algorithmic design can be found in this [manuscript](https://www.biorxiv.org/content/10.1101/2022.08.02.502537v1).

EQO mainly serves as a "hypothesis generating tool" for microbiome studies. By pinpointing a group of species with strongest statistical association with a functional readout (substrate consumption, metabolite production, pathogen inhibition, etc.), EQO can inform further mechanistic studies regarding physiological properties and functional roles of the predicted species. This can be in particular useful for systems where annotation-based approaches are limited.

## 2. Installation

*mEQO* can be installed by 

```R
devtools::install_github("Xiaoyu2425/mEQO")
```

These R packages will be automatically imported as dependencies of mEQO: GA, doParallel, data.table, reshape2. In addition, to implement EQO with the boolean least square (BLS) algorithm, a Gurobi optimizer and its R interface will be required as dependencies. A free license for academic use can be easily acquired on their [website](https://www.gurobi.com). An installation guide can be found [here](https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation_guide.html). 

## 3. A quick example 

Input for EQO is simply the microbiome compositon and the functional/environmental readout of interest. Microbiome composition should be a matrix whose rows are samples and columns are taxa. Column names should be provided. Functional/environmental variable should be a vector whose length is the number of samples. You can find an example (`Microbiome` and `trait`) for input format after loading mEQO by running `library(mEQO)`.

Here we demonstrate the use of mEQO with this example of a synthetic dataset. As shown in the figure below, in this mock microbiome, the relative abundance of single species (color strips) is uncorrelated with the functional trait (black dots). The goal of the algorithm is to identify a group of species whose total abundance is most strongly correlated with the functional trait. 

<img src="https://tva1.sinaimg.cn/large/008vxvgGgy1h7dniv5szcj31900u0whr.jpg" alt="composition1" height="300" width="450" />

To that end, we can run one of the functions below to implement EQO.  

```R
EQO_bls(Microbiome,trait) # implemented by BLS
EQO_ga("c",Microbiome,trait,maxIter=100) # or, implemented by GA
```

We have provided two different algorithms for implementing EQO, i.e., Boolean Least Square (BLS) and Genetic Algorithm (GA), each with its advantages and disadvantages (detailed in the next section). In the above example of GA, "c" indicates a continuous trait variable. In addition, `EQO_ga()` can also support a trait variable that is discrete or uniform. You can specify "d" for a discrete trait variable (e.g., healthy or diseased hosts) or "u" for a uniform trait variable (e.g., stable community composition across samples).

Here, the best group of species found by the algorithm that is coupled to the trait variable is comprised of species 1, 2 and 3, as shown in the figure below. 

<img src="https://tva1.sinaimg.cn/large/008vxvgGgy1h7dniyrywvj31900u0dhz.jpg" alt="composition2" height="300" width="450" />

For microbiome studies, in addition to obtaining a group of species with strongest statistical power, we may be more interested in understanding which single species or what combinations of species within that group have the highest importance. 

We have two suggestions in this regard. 

Firstly, you can apply regularization to the algorithm by specifying *Nmax* as the maximal number of species allowed in the group. In statistics, this is often based on Akaike Information Criterion (AIC). Here, the minimal AIC is achieved at *Nmax=3*, as shown in the figure below. 

```R
aic<-sapply(1:5,function(N){
	assemblage<-EQO_bls(Microbiome,trait,Nmax=N)$abundance
	return(AIC(lm(trait~assemblage))+2*(N-1))
})
```

<img src="https://tva1.sinaimg.cn/large/008vxvgGgy1h7dexjolgcj31900u0wfy.jpg" alt="aic" height="280" width="400" align="center" />

Secondly, we recommend illustrating your functional group with a Cross-validation-based Aggregation Network (CAN). We provide a `CAN()` function in the *mEQO* package for that. In short, this function performs the following: (1) splitting your original dataset into a test set and a validation set, (2) performing EQO with the test set and cross-validation with the validation set, 3) evaluating relative importance of individual species or pairs by repeating (1) and (2) for multiple times. Remarkably, CAN can illustrate patterns that are very different from a traditional association network such as a co-occurrence network. In CAN, each nodes still denotes a species, while the weight of an edge indicates the strength of cross-validated correlation with the functional trait when the two connected species are coarse-grained in a group. In other words, species nodes connected by a stronger edge are more likely to co-exist in a functionally cohesive group. As is shown in the following figure, species 1, 2, and 3 stand out from the cross-validation as a highly interconnected module, strongly suggesting their emergent ecological role as a group together. 

```R
can<-CAN("ga_c",Microbiome,trait,maxIter=100,tm=10)

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

Finally, `EQO_ga()` also supports cases where you want certain species to be included in the functional group, when you have *a priori* knowledge of the functional role of that species. Then mEQO can start on the basis of that partially known functional group and continue optimizing by combining other species. In the package, this can be implemented by specifying the argument *pk* for function `EQO_ga()`. As a vector whose length is equal to the total number of species in the microbiome, *pk* indicates partially known functional group based on a priori knowledge, with 1 for species forced to be included in the targeted group and 0 for the other unknown species. EQO will then search on the basis of the provided partially known group without removing the designated species. For instance, in the case where we want to ensure species 4 is included in the final group, we would include the following argument:

```R
EQO_ga("c",Microbiome,trait,pk=c(0,0,0,1,0,0,0,0),maxIter=100)
```

## 4. Important notes on algorithms

Coarse-graining of species in a microbiome, including EQO, is a combinatorial problem in nature. Unfortunately, combinatorial optimization is still one of the major open challenges in modern operational research, especially for large-scale problems. Currently, there are three possible approaches to tackle the combinatorial optimization problem in EQO, with different advantages and disadvantages.

i) EQO can be reformulated as a Boolean Least Square (BLS) problem when the functional trait is continuous. We can leverage a mature commercial mixed interger programming solver, the Gurobi optimizer, to solve the BLS problem accurately. Due to its combinatorial nature, it is crucial to reduce the search space in order to solve BLS efficiently (There are in total 75,287,520 possible ways to pick 5 species from 100 species). This can be done by either decreasing the total number of species in the microbiome `n` by filtering rare species prensent only in very few samples, or limiting the maximal number of species allowed to be included in the group by specifiying `Nmax` in functions. The computational performance of BLS can be dependent on specific data structure, but we generally recommend `n<200` and `Nmax<10` for a personal laptop, so that BLS can be solved within minutes.  

ii) Genetic algorithm (GA) can be more suitable for larger-scale problems. [GA](https://cran.r-project.org/web/packages/GA/vignettes/GA.html) has been widely used as a successful strategy for stochastic binary searching with convergence in probability. However, the rate of convergence is dependent on size and structure of spefific datasets. Below is an illustration of the accuracy of GA at different number of iterations with simulated microbiome datasets. 

<img src="https://tva1.sinaimg.cn/large/008vxvgGgy1h7hawyk0hrj30pq0roq4r.jpg" alt="截屏2022-10-24 下午9.59.06" height="400" width="380" />

You can determine the optimal number of iterations for your specific datasets based on your computational power and desirable accuracy. Importantly, instead of focusing only on optimality (e.g., having computing time doubled to optimize from r = 0.80 to r = 0.81),  it is of greater ecological significance to ask what species combinations show the most robust pattern across multiple instances of cross-validations. This can be easily done by constructing a CAN and examining strongly-connected nodes, although it might be challenging to accurately solve the exact group with the very optimal performance, in particular for large-scale problems. 

iii) EQO can be also reformulated into an equivalent mixed integer linear programming (MILP) problem, which can be accurately solved by the Gurobi optimizer for continuous, discrete or uniform functional traits. However, this approach only works efficiently for smaller scale problems with limited search space. This approach is not included in the current version of package, but a script can be found [here](https://github.com/Xiaoyu2425/Ensemble-Quotient-Optimization/blob/main/EQO_Reformulation.R). 

## 5. Citation

Please cite [mEQO](https://www.biorxiv.org/content/10.1101/2022.08.02.502537v1.abstract) as an R package, as well as the following dependencies: [GA](https://cran.r-project.org/web/packages/GA/index.html); [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html); [data.table](https://cran.r-project.org/web/packages/data.table/index.html); [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html).

