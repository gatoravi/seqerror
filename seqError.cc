/*
    Estimate Error model based on sequence context.
    Copyright (C) <2013>  <Avinash Ramu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include<iostream>
#include<string>
#include<sstream>
#include<map>
#include<cstdlib>
#include<cmath>

#define RD_CUTOFF 10 //process pileup of lines with RD > RD_CUTOFF
#define MIN_TRIPLET_COUNT 10 // minimum triplet count for a triplet to be output.
using namespace std;

std::map<std::string, float> triplet_cumul_mismatch_counts;
std::map<std::string, double> triplet_cumul_mismatch_fraction;
std::map<std::string, long> triplet_count;
int g_indel_count = 0;

struct line_object{
  string chr;
  char fwd_ref_base;
  char rev_ref_base;
  long pos;
  long read_depth;
  string read_bases;
};

void print_op()
{
  cout<<"Triplet\tTriplet_cumul_mismatch_count\tTriplet_cumul_mismatch_fraction\tTriplet_count\tNormalized_mismatch_fraction\n";
  std::map<std::string, double>::iterator it;
  for(it = triplet_cumul_mismatch_fraction.begin();
      it != triplet_cumul_mismatch_fraction.end();
      it++) {
    if(triplet_count[it->first] > MIN_TRIPLET_COUNT) {
      cout<<it->first;
      cout<<"\t"<<triplet_cumul_mismatch_counts[it->first];
      cout<<"\t"<<it->second<<"\t"<<triplet_count[it->first];
      cout<<"\t"<<it->second/triplet_count[it->first]<<endl;
      //    cout<<"\t"<<exp(it->second/triplet_count[it->first]); //antilog
    }
  }
}


void print_map_kmer_counts()
{
  cout<<"In print map kmer counts."<<endl;
  std::map<std::string, long>::iterator it;
  for(it = triplet_count.begin();
      it != triplet_count.end();
      it++) {
    cout<<it->first<<"\t"<<it->second<<endl;
  }
}



char complement_DNA(char fwd_base)
{
  char rev_base = 'N';
  switch(fwd_base) {
  case 'A':
    rev_base = 'T';
    break;
  case 'a':
    rev_base = 't';
    break;

  case 'C':
    rev_base = 'G';
    break;
  case 'c':
    rev_base = 'g';
    break;

  case 'G':
    rev_base = 'C';
    break;
  case 'g':
    rev_base = 'c';
    break;

  case 'T':
    rev_base = 'A';
    break;
  case 't':
    rev_base = 'a';
    break;

  case 'N':
    rev_base = 'N';
    break;
  default:
    cerr<<"Invalid base: "<<fwd_base;
    exit(EXIT_FAILURE);
  }
  return rev_base;
}



void parse_line_object(string line, line_object& lo1)
{
  stringstream iss(line);
  iss >> lo1.chr;  
  iss >> lo1.pos;

  iss >> lo1.fwd_ref_base;
  if(lo1.fwd_ref_base >= 97)
    lo1.fwd_ref_base -= 32;
#ifdef DEBUG
  cout<<"ref"<<lo1.fwd_ref_base<<endl;
#endif
  iss >> lo1.read_depth;
  iss >> lo1.read_bases;
  lo1.rev_ref_base = complement_DNA(lo1.fwd_ref_base);
}


void process_read_bases(string read_bases, const char* fwd_triplet, char* rev_triplet)
{

  int fwd_match = 0, fwd_mismatch = 0;
  int rev_match = 0, rev_mismatch = 0;
  
#ifdef DEBUG
  cout<<endl<<"Triplet is "<<fwd_triplet<<"\t"<<rev_triplet<<endl;
#endif
  
  if(read_bases.find('+') != std::string::npos || read_bases.find('-') != std::string::npos) {
#ifdef DEBUG
    cout<<endl<<"INDEL "<<read_bases<<" indel_count "<<g_indel_count+1;
#endif
    ++g_indel_count;
    return;
  }

  for(int pos = 0; pos<read_bases.length(); pos++) {
#ifdef DEBUG
    cout<<read_bases[pos];
#endif
    if(read_bases[pos] == '^')
      pos += 1;
    else if(read_bases[pos] == '.')
      fwd_match += 1;
    else if(read_bases[pos] == ',')
      rev_match += 1;
    else if(read_bases[pos] == 'A' || read_bases[pos] == 'C' || read_bases[pos] == 'G' || read_bases[pos] == 'T' )
      fwd_mismatch += 1;
    else if(read_bases[pos] == 'a' || read_bases[pos] == 'c' || read_bases[pos] == 'g' || read_bases[pos] == 't' )
      rev_mismatch += 1;
  }

#ifdef DEBUG
  cout<<endl<<"fwd match "<<fwd_match<<"\t"<<" fwd_mismatch "<<fwd_mismatch;
  cout<<endl<<"rev match "<<rev_match<<"\t"<<" rev_mismatch "<<rev_mismatch<<endl;
#endif

  if(fwd_match + fwd_mismatch != 0 && fwd_mismatch != 0)
     triplet_cumul_mismatch_fraction[fwd_triplet] += float(fwd_mismatch)/float(fwd_match + fwd_mismatch);
    // triplet_cumul_mismatch_fraction[fwd_triplet] += log(float(fwd_mismatch)) - log(float(fwd_match + fwd_mismatch)); /*log space*/
  else
    triplet_cumul_mismatch_fraction[fwd_triplet] += 0;
    // triplet_cumul_mismatch_fraction[fwd_triplet] += -255; /*log space*/
  

  if(rev_match + rev_mismatch != 0 && rev_mismatch !=0)
    triplet_cumul_mismatch_fraction[rev_triplet] += float(rev_mismatch)/float(rev_match + rev_mismatch);
    // triplet_cumul_mismatch_fraction[rev_triplet] += log(float(rev_mismatch)) - log(float(rev_match + rev_mismatch)); /*log space*/
  else 
    triplet_cumul_mismatch_fraction[rev_triplet] += 0;
    // triplet_cumul_mismatch_fraction[rev_triplet] += -255; /*log space*/

  triplet_count[fwd_triplet] += 1;
  triplet_count[rev_triplet] += 1;

  triplet_cumul_mismatch_counts[fwd_triplet] += fwd_mismatch;
  triplet_cumul_mismatch_counts[rev_triplet] += rev_mismatch;
  

}



void process_mpileup()
{
  string current_line, next1_line, next2_line;
  line_object current_line_object, next1_line_object, next2_line_object;
  char  fwd_triplet[3], rev_triplet[3];
  string old_chr = "NA";
  long first_pos = 0, last_pos = 0;
  int line_number = 0;
  bool first_line = true;

  if(!getline(cin, current_line))
    exit(EXIT_FAILURE);
  if(!getline(cin, next1_line)) 
    exit(EXIT_FAILURE);


  while(getline(cin, next2_line)) {

    line_number++;
    
    parse_line_object(current_line, current_line_object);
    if(first_line) {
      first_pos = current_line_object.pos;
      first_line = false;
    }
    parse_line_object(next1_line, next1_line_object);
    parse_line_object(next2_line, next2_line_object);
#ifdef DEBUG    
    cout<<current_line_object.fwd_ref_base<<endl;
#endif
    if(current_line_object.chr != old_chr)
      line_number = 1;
    if(line_number > 2) {
      fwd_triplet[2]= current_line_object.fwd_ref_base;
      rev_triplet[2]= current_line_object.rev_ref_base;
      rev_triplet[1]= next1_line_object.rev_ref_base;
      rev_triplet[0]= next2_line_object.rev_ref_base;
      fwd_triplet[3] = '\0';
      rev_triplet[3] = '\0';
      if(current_line_object.read_depth > RD_CUTOFF)
	process_read_bases(current_line_object.read_bases, fwd_triplet, rev_triplet);
      fwd_triplet[0] = fwd_triplet[1];
      fwd_triplet[1] = fwd_triplet[2];
    }
    else {
      fwd_triplet[line_number - 1] = current_line_object.fwd_ref_base;
    }
    old_chr = current_line_object.chr;
    current_line = next1_line;
    next1_line = next2_line;
  };
  last_pos = current_line_object.pos;
  cerr<<"first_pos: "<<first_pos<<"\tlast_pos: "<<last_pos<<endl;
}




int main(int argc, char* argv[])
{
  process_mpileup();
  print_op();
  //  print_map_kmer_counts();
  cerr<<endl<<"The total number of indels is "<<g_indel_count;
  cerr<<endl<<"Done!";
  exit(EXIT_SUCCESS);
}
