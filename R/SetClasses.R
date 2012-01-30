#TODO:to add comments for each slot
setClass("peptideSet", contains=c("ExpressionSet")
					, representation(
#									featureSequence="character"#better to merge to assay data?
#									, featureID="character"#better to  merge to assay data?
#									, featureAnnotation="data.frame" #antibody binding information
#									,featurePosition="integer"#AA position in HIV envolope
									featureRange="RangedData"
									)
					,prototype(
#								featureSequence=character(0)
#								,featureID=character(0)
#								, featureAnnotation=data.frame()
#								,featurePosition=integer(0)
#								featureRange=RangedData()
								)
#	,validity=function(object){		
#		if(!all(is.character(object@featureSequence)))
#		{
#		stop ("All featureSequence must be a character")	
#		}
#				
#		if(!is.data.frame(object@featureAnnotation))
#		{
#		stop ("featureAnnotation must be a data.frame")	
#		}		
#
#		if(!all(is.integer(object@featurePosition)))
#		{
#		stop ("All featurePosition must be a integer")	
#		}		
#	}
	)

