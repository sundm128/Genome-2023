metabolites <- read.delim('metabolites.txt', row.names = 1)
group <- read.delim('group.txt', row.names = 1)
metabolites <- data.frame(t(metabolites))
metabolites <- metabolites[rownames(group), ]
metabolites <- cbind(metabolites, group)
set.seed(123)
train <- sample(nrow(metabolites), nrow(metabolites)*0.7)
metabolites_train <- metabolites[train, ]
metabolites_test <- metabolites[-train, ]
library(randomForest)
set.seed(123)
metabolites_train.forest <- randomForest(group~., data = metabolites_train, importance = TRUE)
metabolites_train.forest
train_predict <- predict(metabolites_train.forest, metabolites_train)
compare_train <- table(train_predict, metabolites_train$group)
compare_train
sum(diag(compare_train)/sum(compare_train))
test_predict <- predict(metabolites_train.forest, metabolites_test)
compare_test <- table(metabolites_test$group, test_predict, dnn = c('Actual', 'Predicted'))
compare_test
importance_metabolites <- data.frame(importance(metabolites_train.forest))
head(importance_metabolites)
importance_metabolites <- importance_metabolites[order(importance_metabolites$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_metabolites)
write.table(importance_metabolites, 'importance_metabolites.txt', sep = '\t', col.names = NA, quote = FALSE)
varImpPlot(metabolites_train.forest, n.var = min(20, nrow(metabolites_train.forest$importance)), main = 'Top 20 - variable importance')
set.seed(123)
metabolites_train.cv <- replicate(5, rfcv(metabolites_train[-ncol(metabolites_train)], metabolites_train$group, cv.fold = 10,step = 1.5), simplify = FALSE)
metabolites_train.cv
metabolites_train.cv <- data.frame(sapply(metabolites_train.cv, '[[', 'error.cv'))
metabolites_train.cv$otus <- rownames(metabolites_train.cv)
metabolites_train.cv <- reshape2::melt(metabolites_train.cv, id = 'metabolite')
metabolites_train.cv$metabolites <- as.numeric(as.character(metabolites_train.cv$otus))
library(ggplot2)
library(splines)
p <- ggplot(metabolites_train.cv, aes(metabolites, value)) +
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of metabolites', y = 'Cross-validation error')
p

