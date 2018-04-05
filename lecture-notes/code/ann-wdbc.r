# load our data set
wdbc <- read.csv(file="wdbc.csv", header=TRUE, sep=",")

# make sure the ANN library is available
library(neuralnet)

# convert the labels into numeric labels and put them into a data frame
Diagnosis.numeric <- as.numeric(wdbc$Diagnosis)
wdbc.df <- data.frame(wdbc,Diagnosis.numeric)

# train a tree
net <- neuralnet(
    Diagnosis.numeric ~ radius1+texture1+perimeter1+area1+smoothness1+compactness1+concavity1+concave_points1+symmetry1+fractal_dimension1,
    wdbc.df,
    threshold=0.01,
    stepmax="1000000",
    lifesign="none",
    hidden=5)

# display the ANN
plot(net)

# the training predictions from the ANN are numeric values,
# turn them into labels by rounding
predicted.labels <- round(net$net.result[[1]])

# plot the confusion matrix
print(table(iris.df$Species.numeric,predicted.labels))

