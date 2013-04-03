#TODO:to add comments for each slot
setClass("peptideSet",
    contains=c("ExpressionSet"),
    representation(featureRange="RangedData")
)

