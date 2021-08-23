make_all_seqs <- function(nts, l){
  seqs <- nts
  for(i in 1:(l-1)){
    new_seqs <- c()
    for(seq in seqs){
      for(nt in nts){
        new_seqs <- c(new_seqs, paste0(seq, nt))
      }
    }
    seqs <- new_seqs
  }
  return(new_seqs)
}

check_orthogonal <- function(seq, bcds, min_different){
  lowest_mm <- 99999999 # mismatches - assume orthogonal to start
  for(bcd in bcds){
    mm <- 0
    for(j in 1:str_length(seq)){
      if(!str_sub(bcd,j,j)==str_sub(seq,j,j)){
        mm <- mm+1
      }
    }
    if(mm < lowest_mm){
      lowest_mm <- mm
    }
  }
  if(lowest_mm >= min_different){
    return(T)
  } else {
    return(F)
    
  }
}

gen_barcodes <- function(n_bcds,
                         bc_l,
                         min_different){
  nts <- c("A", "C", "G", "T")
  final_bcds <- c()
  
  all_seqs <- sample(make_all_seqs(nts, bc_l))
  
  bcds <- c(all_seqs[1])
  
  for(i in 2:length(all_seqs)){
    if(check_orthogonal(all_seqs[i], bcds, min_different)){
      bcds <- c(bcds, all_seqs[i])
      if(length(bcds) == n_bcds){
        break
      }
    }
  }
  
  return(bcds)
  
}

gen_barcodes(4, 4, 4)
