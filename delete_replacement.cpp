#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include  "omp.h"
#include "stdio.h"

using namespace std;
using namespace Eigen;

char nucleotide_table[4] = { 'A', 'C', 'G', 'T' }; // A C G T;

struct node {
	vector<string> xRcx;
	char row_index;
};

struct rnode {
	vector<string> set;
	int flag;
};

// Decimal integer to binary character string.
string toBinary(long long  value, int words_length) {
	string words;
	while (value != 0) {
		words += (value % 2 == 0 ? '0' : '1');
		value = value / 2;
	}
	reverse(words.begin(), words.end());
	if (words.size() < words_length)
		words.insert(0, words_length - words.size(), '0');
	return words;
}

// Get the reverse-complement of x. 
string get_RC_x(string x) {
	for (int i = 0; i < x.size(); i++) {
		switch (x[i]) {
		case 'A': x[i] = 'T'; break;
		case 'C': x[i] = 'G'; break;
		case 'G': x[i] = 'C'; break;
		case 'T': x[i] = 'A'; break;
		default: cout << "Invalid character";
		}
	}
	reverse(x.begin(), x.end());
	return x;
}

// Recursion produces a permutation of the specified length.
void dfs(int times, int m, unordered_set<string>& findRepeat, vector<node>& m_subsequences, string& s) {

	for (int i = 0; i < 4; i++) {
		s[times] = nucleotide_table[i];
		if (times == m - 1) {
			if (findRepeat.find(s) == findRepeat.end()) {
				struct node sub_node;
				if (s != get_RC_x(s)) {
					sub_node.xRcx = { s , get_RC_x(s) };
					m_subsequences.push_back(sub_node);
					findRepeat.emplace(s);
					findRepeat.emplace(get_RC_x(s));
				}
			}
		}
		else dfs(times + 1, m, findRepeat, m_subsequences, s);
	}
}

// Generate x, RC_x sequence pair.
void get_Permutation(vector<node>& m_subsequences, int m) {
	vector<string> x, RC_x;
	unordered_set<string> findRepeat;
	string s(m, '0');
	dfs(0, m, findRepeat, m_subsequences, s);
}

// Output all xRcx.
void test_output(vector<node>& m_subsequences) {
	// output DNA m_subsequences.
	cout << "m_subsequences��" << endl;
	for (int i = 0; i < m_subsequences.size(); i++)
		cout << m_subsequences[i].xRcx[0] << " ";
	cout << endl;
	for (int i = 0; i < m_subsequences.size(); i++)
		cout << m_subsequences[i].xRcx[1] << " ";
	cout << endl;

}

// Get all replace of string.
string get_replace(string x, string rule) {
	vector<string> set;
	for (int i = 0; i < x.size(); i++) {
		switch (x[i]) {
		case 'A': x[i] = rule[0]; break;
		case 'T': x[i] = rule[1]; break;
		case 'C': x[i] = rule[2]; break;
		case 'G': x[i] = rule[3]; break;
		default: cout << "Invalid character";
		}
	}
	string s1 = "";
	for (int i = 0; i < x.size(); i++) {
		if (s1.size() == 2) {
			set.push_back(s1);
			s1 = "";
		}
		s1 += x[i];
	}
	set.push_back(s1);
	sort(set.begin(), set.end());
	string s2;
	for (int i = 0; i < set.size(); i++)
		s2 += set[i];
	return s2;
}

// Delete all replace sequence equal to some one of ALl_Set.
void replace_function(unordered_map<string,int>& allset) {
	vector<string> rules = {"TACG","TAGC","ATGC","CGAT","CGTA","GCAT","GCTA"};
	int number2 = rules.size();
	int  flagg = 0;
	for (auto i = allset.begin(); i != allset.end(); i++) {
		if (i->second == 0) continue;
		for (int j = 0; j < number2; j++) {
			string s = get_replace(i->first, rules[j]);
			//string s = "";
			if (allset.find(s) != allset.end()) {
				allset[s] = 0;
				flagg++;
			}
		}
	}
	cout << "-------------------------------------------------------------------" << endl;
	for (auto i = allset.begin(); i != allset.end(); i++) {
		if (i->second == 1) {
			int times = 0;
			cout << "{";
			for (int j = 0; j < i->first.size(); j++) {
				if (j  % 2 == 0 && j != 0) cout << ",";
				cout << i->first[j];
			}
			cout << "}";
			cout << endl;
		}	
	}

	cout << flagg;
}

// Get the set of all possible combinations for x and Rc(x).
void all_set(vector<node>& m_subsequences) {
	int number_order = m_subsequences.size();
	long long flag = pow(2, number_order) - 1;
	unordered_map<string, int> allset;
	vector<string> set;
	int number = 0;
	for (long long times = 0; times <= flag; times++) {
		string s;
		string free_RC_set = toBinary(times, m_subsequences.size());
		for (int i = 0; i < free_RC_set.size(); i++) {
			set.push_back(m_subsequences[i].xRcx[free_RC_set[i] - '0']);
		}
		number++;
		sort(set.begin(), set.end());
		for (int i = 0; i < set.size(); i++)
			s += set[i];
		allset.emplace(s, 1); 
		s = ""; set = {};
	}
	/*for (auto i = allset.begin(); i != allset.end(); i++) {
			cout << i->first << i->second <<  endl;
	}*/
	//cout << number;
	replace_function(allset);
}


int main() {

	int m = 2;
	long long times, flag;
	string prefix_string;
	vector<node> m_subsequences;
	vector<node> binary_m_subsequences;
	get_Permutation(m_subsequences, m);
	test_output(m_subsequences);
	all_set(m_subsequences);

	return 0;
}
