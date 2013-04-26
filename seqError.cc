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

using namespace std;

std::map<std::string, float mismatch_fraction> triplet_cumul_mismatch_fraction;
//std::map<std::string, int count> triplet_count;

void process_read_bases(string read_bases, const char* triplet)
{
  cout<<endl<<"Triplet is "<<triplet;
  for(int pos = 0; pos<read_bases.length(); pos++) {
    cout<<read_bases[pos];
  }
}


void process_mpileup()
{
  string input_line;
  char  triplet[3];
  string chr, read_bases;
  char ref_base;
  long pos, read_depth;
  int line_number = 0;
  while(getline(cin, input_line)) {
    line_number++;
    stringstream iss(input_line);
    iss >> chr;
    iss >> pos;
    iss >> ref_base;
    iss >> read_depth;
    iss >> read_bases;
    cout<<ref_base<<endl;
    if(line_number > 2) {
      triplet[2]= ref_base;
      triplet[3] = '\0';
      std::cout << "inside avi's program read_bases: "<< read_bases  << endl;
      process_read_bases(read_bases, triplet);
      cin.clear();
      triplet[0] = triplet[1];
      triplet[1] = triplet[2];
    }
    else
      triplet[line_number - 1] = ref_base;
  };
}

int main(int argc, char* argv[])
{
  process_mpileup();
}
