/*
 * tfm_index_invert.cpp for Edge minimization in de Bruijn graphs
 * Copyright (c) 2019 Uwe Baier, Pascal Weber All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <iostream>
#include <map>
#include <deque>
#include <utility>
#include <vector>
#include <cstdio>

#include "tfm_index.hpp"

using namespace std;
using namespace sdsl;

typedef typename sdsl::int_vector<>::size_type size_type;

void printUsage( char **argv ) {
	cerr << "USAGE: " << argv[0] << " TFMFILE" << endl;
	cerr << "TFMFILE:" << endl;
	cerr << "  File where to store the serialized trie" << endl;
};

int main(int argc, char **argv) {
	//check parameters
	if (argc < 2) {
		printUsage( argv );
		cerr << "At least 1 parameter expected" << endl;
		return 1;
	}
	
	//load tunneled fm index
	tfm_index<> tfm;
        construct_from_pfwg(tfm, argv[1]);

	char *S = new char[tfm.size()];
	auto p = tfm.end();
	for (size_type i = 0; i < tfm.size(); i++) {
            char c = (char)tfm.backwardstep( p );
            S[tfm.size() - i - 1] = c;
            //cout << "\t pos:" << p.first << " c: " << c << endl;
	}
        /*
        cout << "\t";
	for (size_type i = 0; i < tfm.size(); i++) {
            cout << S[i];
        }
        cout << endl;
        */
        
        char *name;
        int e = asprintf(&name, "%s.%s", argv[1], "untunneled");
        if (e == -1) {cout << "ERROR while creating the output file name!" << endl; return 1; }
        FILE *fout = fopen(name, "w");
        fwrite(S, sizeof(char), tfm.size(), fout);
        fclose(fout);
}
