# load our data set
data(iris)

# make sure the rpart library is available
library(rpart)

# train a tree
tree <- rpart(
    Species ~ .,
    data=iris,
    maxdepth=10,
    minsplit=1,
    minbucket=1,
    cp = 0.0001,
    method="class")


# display the tree
plot(tree)
text(tree)

# confusion matrix
model.output = predict(tree, type="class")
print(table(iris$Species, model.output))


