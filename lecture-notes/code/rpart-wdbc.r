# load our data set
wdbc <- read.csv(file="wdbc.csv", header=TRUE, sep=",")

# make sure the rpart library is available
library(rpart)

# train a tree
tree <- rpart(
    Diagnosis ~ radius1+texture1+perimeter1+area1+smoothness1+compactness1+concavity1+concave_points1+symmetry1+fractal_dimension1,
#    Diagnosis ~ .,
    data=wdbc,
    maxdepth=30,
    minsplit=1,
    minbucket=1,
    cp = 0.01,
    method="class")


# display the tree
plot(tree)
text(tree)

# confusion matrix
model.output = predict(tree, type="class")
print(table(wdbc$Diagnosis, model.output))

