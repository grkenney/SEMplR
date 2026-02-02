# reverse complement SEMs
reverseComplementSEM <- function(sem) {
  if(!is(sem, "SNPEffectMatrix")) {
    rlang::abort("sem must be an object of class SNPEffectMatrix")
  }
  sem_matrix <- getSEM(sem)
  
  # reverse the SEM
  r_sem_matrix <- sem_matrix[seq(nrow(sem_matrix), 1, -1)]
  
  # complement the SEM
  # make sure the columns are ordered alphabetically
  rc_sem_matrix <- r_sem_matrix[, c("A", "C", "G", "T")]
  # do the complement
  colnames(rc_sem_matrix) <- c("T", "G", "C", "A")
  # reorder the columns alphabetically
  rc_sem_matrix <- rc_sem_matrix[, c("A", "C", "G", "T")]
  
  rc_sem <- SNPEffectMatrix(rc_sem_matrix, 
                            getBaseline(sem), 
                            semId = getSEMId(sem))
  
  return(rc_sem)
}