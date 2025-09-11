#include <Rcpp.h>
using namespace Rcpp;

float score_seq(Rcpp::NumericMatrix sem,
                     Rcpp::StringVector dna_sequence) {
  Rcpp::CharacterVector sem_cols = colnames(sem);
  // Rcpp::Rcout << sem_cols << "\n";
  
  float score = 0.00;
  // iterate over rows of the SEM
  for(int i = 0; i<sem.nrow(); i++){
    char base = dna_sequence[0][i];
    // iterate over columns to find matching base score
    for(int j = 0; j<sem.ncol(); j++)
      if (base == sem_cols[j][0]){
        score += sem(i, j);
      }
  }
  return score;
}

//[[Rcpp::export]]
Rcpp::List score_frames(Rcpp::NumericMatrix sem,
                        Rcpp::StringVector dna_sequence,
                        int first_frame_index,
                        int last_frame_index) {
  int sem_nrows = sem.nrow();
  int n_frames = last_frame_index - first_frame_index + 1;
  Rcpp::IntegerVector frame_starts = Rcpp::seq(first_frame_index, 
                                               last_frame_index);
  
  float neg_inf = -std::numeric_limits<float>::infinity();
  Rcpp::List max_score = Rcpp::List::create(Rcpp::Named("score") = neg_inf,
                                            Rcpp::_["index"] = 0,
                                            Rcpp::_["seq"] = "");
  
  
  for(int i = 0; i<n_frames; i++) {
    std::string cppString = Rcpp::as<std::string>(dna_sequence);
    Rcpp::StringVector frame = cppString.substr(frame_starts[i], sem_nrows);
    float frame_score = score_seq(sem, frame);
    
    // Rcpp::Rcout << "frame: " << frame << "\n";
    // Rcpp::Rcout << "score: " << frame_score << "\n";
    
    float currentMax = Rcpp::as<float>(max_score["score"]); 
    if(frame_score > currentMax) {
      max_score["score"] = frame_score;
      max_score["index"] = frame_starts[i];
      max_score["seq"] = frame;
    }
  }
  return max_score;
}


//[[Rcpp::export]]
Rcpp::DataFrame scoreSequence(Rcpp::NumericMatrix sem, 
                           Rcpp::CharacterVector dna_sequences,
                           int nFlank,
                           float bl,
                           Rcpp::CharacterVector seqIds) {
  int n_seqs = dna_sequences.size();
  
  Rcpp::NumericVector score(n_seqs);
  Rcpp::NumericVector scoreNorm(n_seqs);
  Rcpp::IntegerVector index(n_seqs);
  Rcpp::CharacterVector seq(n_seqs);
  Rcpp::CharacterVector seqId(n_seqs);

  Rcpp::List max_score = Rcpp::List::create(Rcpp::Named("score") = 0,
                                            Rcpp::_["scoreNorm"] = 0,
                                            Rcpp::_["index"] = 0,
                                            Rcpp::_["seq"] = "",
                                            Rcpp::_["seqId"] = "");
  
  for (int i = 0; i<n_seqs; i++) {
    std::string dna_sequence = Rcpp::as<std::string>(dna_sequences[i]); 
    int first_frame_index = nFlank + 1 - sem.nrow();
    
    int last_frame_index = 0;
    // if no flank, score up to last possible index given sem length
    // if flank, score up to last index including one non-flank character
    if (nFlank == 0) {
      last_frame_index = dna_sequence.length() - sem.nrow();
    } else {
      last_frame_index = dna_sequence.length() - nFlank - 1;
    }
    
    //Rcpp::Rcout << "nFlank: " << nFlank << "\n";
    //Rcpp::Rcout << "length: " << dna_sequence << "\n";
    //Rcpp::Rcout << "length: " << dna_sequence.length() << "\n";
    //Rcpp::Rcout << "last frame: " << last_frame_index << "\n";
    
    
    // if the number of characters from last index is less than the sem length
    // move index back to last possible index given sem length
    if (dna_sequence.length() - last_frame_index < sem.nrow()) {
      //Rcpp::Rcout << "adjusting last frame\n";
      last_frame_index = dna_sequence.length() - sem.nrow();
      //Rcpp::Rcout << "last frame: " << last_frame_index << "\n";
    }
    
    if (first_frame_index < 0) {
      first_frame_index = 0;
    }
    
    //Rcpp::Rcout << "first frame: " << first_frame_index << "\n";
    //Rcpp::Rcout << "last frame: " << last_frame_index << "\n";
    
    max_score = score_frames(sem, dna_sequence, 
                             first_frame_index, 
                             last_frame_index);
    score[i] = max_score["score"];
    scoreNorm[i] = (std::pow(2, score[i]) - std::pow(2, bl)) / std::pow(2, bl);
    index[i] = int(max_score["index"]) + 1;
    std::string max_seq = Rcpp::as<std::string>(max_score["seq"]); 
    seq[i] = max_seq;
    seqId[i] = Rcpp::as<std::string>(seqIds[i]);
    }
  
  Rcpp::DataFrame fss = Rcpp::DataFrame::create(
    Rcpp::Named("score") = score,
    Rcpp::Named("scoreNorm") = scoreNorm,
    Rcpp::Named("index") = index,
    Rcpp::Named("seq") = seq,
    Rcpp::Named("seqId") = seqId
  );
  
  return fss;
}
