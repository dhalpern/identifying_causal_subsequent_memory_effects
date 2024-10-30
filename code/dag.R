library(dagitty)
library(rethinking)
dag_concreteness <- dagitty("dag { 
                   C -> CNA -> M 
                   C -> NCNA
                   fMRI <- NCNA
                   }") 
plot(graphLayout(dag_concreteness))
paths(dag_concreteness, "fMRI", "M")
paths(dag_concreteness, "fMRI", "M", c("C"))
adjustmentSets(dag_concreteness, exposure="C", outcome="M")
impliedConditionalIndependencies(dag_concreteness)

dag_fMRI <- dagitty("dag { CNA[unobserved]
                    NCNA[unobserved]
                    EF -> CNA -> M 
                    EF -> NCNA
                    CNA -> fMRI
                    NCNA -> fMRI}") 
plot(graphLayout(dag_fMRI))
paths(dag_fMRI, "fMRI", "M")
paths(dag_fMRI, "fMRI", "M", c("EF"))
impliedConditionalIndependencies(dag_fMRI)
adjustmentSets(dag_fMRI, exposure="fMRI", outcome="M")
drawdag(dag_fMRI)

dag_fMRI2 <- dagitty("dag { CNA[unobserved]
                    NCNA[unobserved]
                    EF -> CNA -> M 
                    EF -> NCNA
                    NCNA -> fMRI}") 
impliedConditionalIndependencies(dag_fMRI2)
adjustmentSets(dag_fMRI2, exposure="fMRI", outcome="M")
drawdag(dag_fMRI2)

dag_fMRI3 <- dagitty("dag { CNA[unobserved]
                    NCNA[unobserved]
                    EF -> CNA -> M 
                    IF -> CNA -> M
                    IF -> NCNA
                    EF -> NCNA
                    NCNA -> fMRI}") 
impliedConditionalIndependencies(dag_fMRI3)
adjustmentSets(dag_fMRI3, exposure="fMRI", outcome="M")
drawdag(dag_fMRI3)
plot(dag_fMRI3)
paths(dag_fMRI3, "fMRI", "M", c("EF"))
