# quick little awk script to process the stdout of the pgp_rcpp() function
# into two different tidy space-delimited files:
#   P_unrel.txt and P_parental.txt
# These are designed to be read in with read_table() and then

BEGIN {
  # Here are the different tag lines that come out:

  # if(WriteRcouts) Rcout << "KPC:"<< k << " " << p << " " << c << " " << IXG(lo, __anck) << " " << IXG(lo, __ancp) << AD(k, 0) << " " << AD(k, 1) << " " << AD(k, 2) << " " << AD(p, 0) << " " << AD(p, 1) << " " << AD(p, 2) << " " << lo << " " << hi << std::endl;
  header["KPC"] = "kIdx pIdx cIdx anc_k anc_p ADk_1 ADk_2 ADk_3 ADp_1 ADp_2 ADp_3 lo hi"

  # if(WriteRcouts) Rcout << "PAk_un: "<< PAk_un << std::endl;
  header["PAk_un"] = "PAk_un"

  # if(WriteRcouts) Rcout << "PGk_un: "<< m << " " << midx << " " << isD(midx)  << " "<< g << " " << gpar << " " << geno_prob << std::endl;
  header["PGk_un"] = "m midx isDiag gk gp geno_prob"

  # if(WriteRcouts) Rcout << "PAk_par: " << PAk_par << std::endl;
  header["PAk_par"] = "PAk_par"
  
  
  # if(WriteRcouts) Rcout << "AsApAn: " << As << " " << An << " " << An << std::endl;
  header["AsApAn"] = "As Ap An"

  #if(WriteRcouts) Rcout << "PGk_par: "<< m << " " << midx << " " << isDiag(midx)  << " " << f1s << " " << f0s << " " << f1p << " " << f0p << " " << gk << " " << gp << " " << geno_prob << std::endl;
  header["PGk_par"] = "m midx isDiag f1s f0s f1p f0p gk gp geno_prob"

  # go ahead and print the headers
  fileUn = "P_unrel.txt"
  filePar = "P_parental.txt"
  
  FS = ":"  # colon is the input field separator for us
  OFS = " "  # space is the output field separator for us

  print header["KPC"], header["PAk_un"], header["PGk_un"] > "P_unrel.txt";
  print header["KPC"], header["PAk_par"], header["AsApAn"], header["PGk_par"] > "P_parental.txt";
}


/KPC:/ {KPC=$2}
/PAk_un:/ {PAk_un=$2}
/PGk_un:/ {print KPC, PAk_un, $2 > "P_unrel.txt"}
/AsApAn:/ {AsApAn = $2}
/PAk_par:/ {PAk_par = $2}
/PGk_par:/ {print KPC, AsApAn, PAk_par, PGk_par, $2 > "P_parental.txt"}



