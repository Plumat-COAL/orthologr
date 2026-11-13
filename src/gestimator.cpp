/*

Copyright (C) 2003-2009 Kevin Thornton, krthornt[]@[]uci.edu

Remove the brackets to email me.

This file is part of libsequence.

libsequence is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libsequence is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
long with libsequence.  If not, see <http://www.gnu.org/licenses/>.


Modified by Sarah Scharfenberg and Hajk-Georg Drost 2014 to work 
in orthologr without using external libraries from libsequence.

All changes are also free under the terms of GNU General Public License
version 2 of the License, or any later version.

*/

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <utility>      // std::pair, std::make_pair
#include <string>       // std::string
#include <iostream>     // std::cout


#include "Comeron95.h"

using namespace Rcpp;
using namespace std;

// @export
// [[Rcpp::export]]
void gestimator(std::string file, std::string file_out = "", int maxHits = 3, 
                bool verbose = false, bool remove_all_gaps = false, bool legacy_saturation = false) {

   if (maxHits < 1 || maxHits > 3)
            Rcpp::stop("maxHits must be 1, 2 or 3.");

   // read codon file -> GetData()
   ifstream in_stream (file.c_str());
   
   // check parameter
   if(!in_stream){
            Rcpp::stop("Input File not found!");
   }
   
   string line;
   vector<pair<string,string> > data;
   string sequence = "";
   pair<string,string> seq;
   bool any_header = false;
   while(getline(in_stream,line)){
           //in_stream >> line;
         
           // files with Windows line ends: the carriage return is not part of the sequence
           while (!line.empty() && (line[line.size() - 1] == '\r' || line[line.size() - 1] == '\n'))
                line.erase(line.size() - 1);
         
           // new seq found: a header line starts with '>' (a '>' further on belongs to the description)
           if (!line.empty() && line[0] == '>'){
                any_header = true;
                // push pevious seq
                if(sequence != "") {
                        seq.second = sequence;
                        data.push_back(seq);
                        sequence = "";
                }
                // get name for next sequence
                // the name is the header up to the first white space (it used to include that space)
                size_t name_end = line.find_first_of(" \t");
                seq.first = (name_end == std::string::npos) ? line.substr(1) : line.substr(1, name_end - 1);
           }
           else {
                // otherwise append to sequence
                sequence.append(line);
           }
   }
   seq.second = sequence;
   data.push_back(seq);
   in_stream.close();

   // the input is checked before anything is computed or written: an alignment that is not a codon alignment
   // is an R error (it used to crash R or to give numbers read from beyond the end of the shorter sequence)
   if (!any_header || data.size() < 2)
            Rcpp::stop("The alignment file must hold at least two sequences in FASTA format.");
   for (size_t i = 0; i < data.size(); i++){
           if (data[i].second.length() != data[0].second.length())
                Rcpp::stop("The aligned sequences '%s' and '%s' have unequal lengths (%d and %d): this is not a codon alignment.",
                           data[0].first, data[i].first, (int) data[0].second.length(), (int) data[i].second.length());
   }
   if (data[0].second.length() < 3)
            Rcpp::stop("The aligned sequences are shorter than one codon (length %d).", (int) data[0].second.length());
   if (data[0].second.length() % 3 != 0)
            Rcpp::stop("The length of the aligned sequences (%d) is not a multiple of 3: this is not a codon alignment.",
                       (int) data[0].second.length());


//  // check correct infile reading 
//   for(int i=0;i<data.size();i++){
//           Rcpp::Rcout << i << endl;
//           Rcpp::Rcout << data[i].first << endl;
//           Rcpp::Rcout << data[i].second << endl;
//           Rcpp::Rcout << data[i].second.length()<<endl;
//   }

   // prepare output
   
   ofstream out_stream;
   if (file_out != "")
        out_stream.open(file_out.c_str());

   if (verbose == true)
    {
      if (file_out != "")
        {
          out_stream <<"seq1\tseq2\tKa\tKs\tomega\tp0\tp2S\tp2V\tp4\tq0\tq2S\tq2V\tq4\t";
          out_stream <<"Aa\tAs\tBa\tBs\tL0\tL2S\tL2V\tL4"<<endl;
        }
      else
        {
          Rcpp::Rcout <<"seq1\tseq2\tKa\tKs\tomega\tp0\tp2S\tp2V\tp4\tq0\tq2S\tq2V\tq4\t";
          Rcpp::Rcout <<"Aa\tAs\tBa\tBs\tL0\tL2S\tL2V\tL4"<<endl;
        }
     }
   
 
   // delete gaps
   if(remove_all_gaps){

		// if (Gapped(&data))
		bool gaps=false;
	    for(size_t i=0;i<data.size();i++){
			if( data[i].second.find("-") != std::string::npos){
				    gaps=true;
			}
	    }
		if(gaps){
		// RemoveGaps(data);

			/*
			  Modifies the data vector to remove all positions that
			  contain the gap character'-'.
			  \param data vector<T> to modify
			*/

		   size_t i, j,datasize=data.size();
		   int length = data[0].second.length ();
		   vector<string> ungapped_sequences(data.size());
		   bool site_is_gapped;

		  for (i = 0; i < length; ++i)
		    {        //iterate over sites
		      for ( j = 0, site_is_gapped = 0;
		            j < datasize;  ++j)
		        {
		          if (data[j].second[i] == '-')
		            {
		              site_is_gapped = 1;
		              j = datasize;
		            }
		        }
		      if (!(site_is_gapped))
		        {
		          for ( j = 0 ; j != data.size();  ++j)
		            ungapped_sequences[j] += data[j].second[i];
		        }
		    }

		  //redo the data
		  for (j = 0; j != datasize ; ++j)
		    {
		      data[j].second = ungapped_sequences[j];
		    }
        
		}
   }
   
   // process ()
   
   RedundancyCom95 code_degeneracy;
   
   // foreach sequence pair in data
   for(size_t i=0;i<data.size();i++){
           for(size_t j=i+1;j<data.size();j++){

                   // an object on the stack: it is freed at the end of every round (it was allocated with new
                   // and never deleted: about 1.5 kB lost per gene pair)
                   Comeron95 comeron(&data[i], &data[j], maxHits, &code_degeneracy, legacy_saturation);
                   Comeron95 *C = &comeron;
                   if (file_out == "")
                        {
                          Rcpp::Rcout.precision(4); 
                          Rcpp::Rcout << data[i].first << '\t';
                          Rcpp::Rcout << data[j].first << '\t';
                          Rcpp::Rcout << C->dn() << '\t';
                          Rcpp::Rcout << C->ds() << '\t';
                          Rcpp::Rcout << C->ratio();
                          if (verbose)
                            {
                              Rcpp::Rcout <<'\t';
                              Rcpp::Rcout <<C->P0() << '\t' <<C->P2S()<<'\t';
                              Rcpp::Rcout <<C->P2V() <<'\t' <<C->P4()<<'\t';
                              Rcpp::Rcout <<C->Q0() << '\t' <<C->Q2S()<<'\t';
                              Rcpp::Rcout <<C->Q2V() << '\t' <<C->Q4()<<'\t';
                              Rcpp::Rcout <<C->aa()<<'\t'<<C->as()<<'\t'<<C->ba()<<'\t'<<C->bs()<<'\t';
                              Rcpp::Rcout.precision(6);
                              Rcpp::Rcout <<C->L0() << '\t' <<C->L2S()<<'\t';
                              Rcpp::Rcout <<C->L2V() << '\t' <<C->L4()<<endl;
                            }
                          else
                            Rcpp::Rcout << endl;
                        }
                      else
                        {
                          out_stream.precision(4);
                          out_stream << data[i].first << '\t';
                          out_stream << data[j].first << '\t';
                          out_stream << C->dn() << '\t';
                          out_stream << C->ds() << '\t';
                          out_stream << C->ratio();
                          if (verbose)
                            {
                              out_stream<<'\t';
                              out_stream <<C->P0() << '\t' <<C->P2S()<<'\t';
                              out_stream <<C->P2V() <<'\t' <<C->P4()<<'\t';
                              out_stream <<C->Q0() << '\t' <<C->Q2S()<<'\t';
                              out_stream <<C->Q2V() << '\t' <<C->Q4()<<'\t';
                              out_stream <<C->aa()<<'\t'<<C->as()<<'\t'<<C->ba()<<'\t'<<C->bs()<<'\t';
                              out_stream.precision(6);
                              out_stream <<C->L0() << '\t' <<C->L2S()<<'\t';
                              out_stream <<C->L2V() << '\t' <<C->L4()<<endl;
                            }
                          else
                            out_stream<<endl;

                         }
                         
                        // manually deleting the previous results
                        //delete C;
                        C=NULL;
                 } // end j
                 
           } // end i
           
           out_stream.close();
} // end fun
