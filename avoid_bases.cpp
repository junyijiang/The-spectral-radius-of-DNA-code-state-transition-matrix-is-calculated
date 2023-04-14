/*
Change the path of the <Eigen/Eigen> package while want running on Linux(my eigen path in linux) .
The second line is the command to run(-O3 let cpp runing with release such that count eigenvalue is faster).
#include </home/environment/eigen/build/usr/local/include/eigen3/Eigen/Eigen>
g++ avoid_bases.cpp -fopenmp -O3 -o a.out && ./a.out
*/
#include <iostream>
#include <unordered_set>
#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include  "omp.h"
#include "stdio.h"

using namespace std;
using namespace Eigen;

char nucleotide_table[4] = { 'A', 'C', 'G', 'T'};

// xRcx contains string x and Rc(x), row_index: 0:select x, 1:select Rc(x), 2: select both of them.
struct node {
	vector<string> xRcx;
	char row_index;
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

// Determines whether the string contains two consecutive and equal numbers of runs(like 0011 or 1100 where 1:A or C, 0:G or T).
bool binary_x(string x) {
	int number = x.size() / 2;
	string run1(number, '0'), run2(number,'1');
	for (int i = 0; i < x.size(); i++) {
		switch (x[i]) {
		case 'A': x[i] = '1'; break;
		case 'C': x[i] = '1'; break;
		case 'G': x[i] = '0'; break;
		case 'T': x[i] = '0'; break;
		default: cout << "Invalid character";
		}
	}
	if (x != run1 + run2 && x != run2 + run1)
		return true;
	else return false;
}

//Judgment whether two sequences are successors.
bool is_Successor(string x_1, string x_2) {
	if (x_1.substr(1) == x_2.substr(0, x_2.size() - 1))
		return true;
	else return false;
}

// Descending order.
bool cmp(double a, double b) {
	return  a > b ? true : false;
}

// Ascending order for node type data.
bool cmp_node(struct node a, struct node b) {
	if (a.row_index == b.row_index)
		return a.xRcx[0] < b.xRcx[0] ? true : false;
	return a.row_index < b.row_index ? true : false;
}

// Judgment whether include a specific base.
bool isbases(string sequences, char base) {
	for (int i = 0; i < sequences.size(); i++) {
		if (sequences[i] == base)
			return false;
	}
	return true;
}

// Recursion produces a permutation of the specified length.
void dfs(int times, int m, unordered_set<string>& findRepeat, vector<node>& m_subsequences, string& s) {

	for (int i = 0; i < 4; i++) {
		s[times] = nucleotide_table[i];
		if (times == m - 1) {
			if (findRepeat.find(s) == findRepeat.end()) {
				struct node sub_node;
				if (s != get_RC_x(s)) {// && binary_x(s) this function delete all sequences like (0011 or 1100).
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

// Generate the successor relation matrix, calculate eigenvalues of the matrix, and output the maximum eigenvalues.
double ES_get_MaxEigenvalue(string& free_RC_set, vector<node>& m_subsequences) {
	int length = free_RC_set.size();
	int index_i, index_j;
	MatrixXd m(length, length);
	//cout << "matrix :" << endl;
	cout << "[";
	for (int i = 0; i < length; i++) {
		index_i = free_RC_set[i] - '0';
		cout << "[";
		for (int j = 0; j < length; j++) {
			index_j = free_RC_set[j] - '0';
			if (is_Successor(m_subsequences[i].xRcx[index_i], m_subsequences[j].xRcx[index_j])) {
				m(i, j) = 1;
				cout << "1";
				if (j != length - 1) cout << ",";
			}
			else {
				m(i, j) = 0;
				cout << "0";
				if (j != length - 1) cout << ",";
			}
		}
		cout << "]";
		if (i != length - 1) {
			cout << ",";
			cout << endl;
		}
	}
	cout << "]" << endl;
	EigenSolver<MatrixXd> es(m);
	VectorXd value = es.eigenvalues().real();
	vector<double> ev;
	for (int i = 0; i < value.size(); i++)
		ev.push_back(value[i]);
	sort(ev.begin(), ev.end(), cmp);
	return ev[0];
}

// AC: 1 GT: 0 Converts all quaternion sequences to binary.
vector<node> fout_to_two(vector<node>& m_subsequences) {
	unordered_set<string> table;
	vector<node> binary_m_subsequences;
	for (int i = 0; i < m_subsequences.size(); i++) {
		string s0, s1;
		struct node binary_node;
		for (int j = 0; j < m_subsequences[i].xRcx[0].size(); j++)
			if (m_subsequences[i].xRcx[0][j] == 'A' || m_subsequences[i].xRcx[0][j] == 'C')
				s0 += '1';
			else s0 += '0';
		for (int j = 0; j < m_subsequences[i].xRcx[1].size(); j++)
			if (m_subsequences[i].xRcx[1][j] == 'A' || m_subsequences[i].xRcx[1][j] == 'C')
				s1 += '1';
			else s1 += '0';
		if (table.find(s0) == table.end() && table.find(s1) == table.end()) {
			table.emplace(s0); table.emplace(s1);
			binary_node.row_index = m_subsequences[i].row_index;
			binary_node.xRcx = { s0, s1 };
			binary_m_subsequences.push_back(binary_node);
		}
	}
	return binary_m_subsequences;
}

// Base distribution rules. when AC quantity equal GT quantity. Select AC in the middle and GT on both sides.
bool m4_isdistribute_AC_TG(string sequences) {
	if (sequences[0] == 'G' || sequences[0] == 'T')
		if (sequences[3] == 'G' || sequences[3] == 'T')
			return  true;
	if (sequences[1] == 'A' || sequences[1] == 'C')
		if (sequences[2] == 'A' || sequences[2] == 'C')
			return  true;
	return false;
}

bool m6_isdistribute_AC_TG(string sequences) {
	if (sequences[0] == 'G' || sequences[0] == 'T')
		if (sequences[5] == 'G' || sequences[5] == 'T')
			return  true;
	return false;
}

int count_bases_number(string sequences, char base1, char base2) {
	int number = 0;
	for (int i = 0; i < sequences.size(); i++)
		if (sequences[i] == base1 || sequences[i] == base2)
			number++;
	return number;
}

// Choose the one with more AC's between x and Rc(x), and there's only one result when m is odd. 
void m_odd_avoid_Bases_AC_TG_balance(vector<node>& m_subsequences) {
	int subset_length = m_subsequences.size();
	for (int i = 0; i < subset_length; i++) {
		if (count_bases_number(m_subsequences[i].xRcx[0], 'A', 'C') > count_bases_number(m_subsequences[i].xRcx[0], 'G', 'T'))
			m_subsequences[i].row_index = '0';
		else if (count_bases_number(m_subsequences[i].xRcx[1], 'A', 'C') > count_bases_number(m_subsequences[i].xRcx[1], 'G', 'T'))
			m_subsequences[i].row_index = '1';
	}
}

// Besides the rule of AC quantity, the rule of base dis tribution is added.
void m4_avoid_Bases_AC_TG_balance(vector<node>& m_subsequences) {
	int subset_length = m_subsequences.size();
	for (int i = 0; i < subset_length; i++) {
		if (count_bases_number(m_subsequences[i].xRcx[0], 'A', 'C') > count_bases_number(m_subsequences[i].xRcx[0], 'G', 'T'))
			m_subsequences[i].row_index = '0';
		else if (count_bases_number(m_subsequences[i].xRcx[1], 'A', 'C') > count_bases_number(m_subsequences[i].xRcx[1], 'G', 'T'))
			m_subsequences[i].row_index = '1';
		else  if (m4_isdistribute_AC_TG(m_subsequences[i].xRcx[0])) {
			m_subsequences[i].row_index = '0';
		}
		else  if (m4_isdistribute_AC_TG(m_subsequences[i].xRcx[1])){
			m_subsequences[i].row_index = '1';			
		}
		else m_subsequences[i].row_index = '2';
	}
	sort(m_subsequences.begin(), m_subsequences.end(), cmp_node);
	for (int i = 0; i < m_subsequences.size(); i++)
		cout << m_subsequences[i].row_index << " ";
	cout << endl;
}

void m6_avoid_Bases_AC_TG_balance(vector<node>& m_subsequences) {
	int subset_length = m_subsequences.size();
	for (int i = 0; i < subset_length; i++) {
		if (count_bases_number(m_subsequences[i].xRcx[0], 'A', 'C') > count_bases_number(m_subsequences[i].xRcx[0], 'G', 'T'))
			m_subsequences[i].row_index = '0';
		else if (count_bases_number(m_subsequences[i].xRcx[1], 'A', 'C') > count_bases_number(m_subsequences[i].xRcx[1], 'G', 'T'))
			m_subsequences[i].row_index = '1';
		else  if (m6_isdistribute_AC_TG(m_subsequences[i].xRcx[0]))
			m_subsequences[i].row_index = '0';
		else  if (m6_isdistribute_AC_TG(m_subsequences[i].xRcx[1]))
			m_subsequences[i].row_index = '1';
		else m_subsequences[i].row_index = '2';
	}
	sort(m_subsequences.begin(), m_subsequences.end(), cmp_node);
	for (int i = 0; i < m_subsequences.size(); i++)
		cout << m_subsequences[i].row_index << " ";
	cout << endl;
}

// Outout all xRcx.
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

// if we have rule to select m_sequence, get selected m_sequences.
string cout_prefix_string(vector<node>& m_subsequences) {
	int subset_length = m_subsequences.size();
	string prefix_string = ""; int prefix_string_length = 0;

	for (int i = 0; i < m_subsequences.size(); i++)
		if (m_subsequences[i].row_index != '2') {
			prefix_string += m_subsequences[i].row_index;
			prefix_string_length++;
		}

	cout << "m_subsequences: " << subset_length << endl：;
	cout << "selected_sequences: " << prefix_string_length << endl：;
	cout << "unseleted_sequences: " << subset_length - prefix_string_length << endl;
	//cout << prefix_string << " " << prefix_string.size() << endl;

	return prefix_string;
}

// test odd m or one sequences
void test_one_sequences(FILE* data, string& prefix_string, vector<node>& m_subsequences) {
	
	string free_RC_set = prefix_string;
	cout << "best_sequence: " << endl;
	for (int i = 0; i < free_RC_set.size(); i++)
		cout << m_subsequences[i].xRcx[free_RC_set[i] - '0'] << " ";
	cout << endl;
	//cout << "state transtion matrix:" << endl;
	double max_eigenvalue = ES_get_MaxEigenvalue(free_RC_set, m_subsequences);
	fprintf(data,"max_eigenvalue: %.18lf \n", max_eigenvalue);
	printf("max_eigenvalue: %.18lf \n", max_eigenvalue);

}

// test even m
void test_all_sequences(FILE* data, long long& times, long long& flag, string prefix_string, vector<node>& m_subsequences) {
	
	int number_order = m_subsequences.size() - prefix_string.size();
	cout << prefix_string.size();
	flag = pow(2, number_order) - 1;
	cout <<  "Cout times: " << flag << endl;
	//set 64 threads and record iteration time.
	#pragma omp parallel for num_threads(64)
	for (times = 0; times <= flag; times++) {
	string free_RC_set = prefix_string + toBinary(times, number_order);
	//string free_RC_set = toBinary(times, m_subsequences.size());
	double max_eigenvalue = ES_get_MaxEigenvalue(free_RC_set, m_subsequences);
		if (max_eigenvalue > 2.246) {
			for (int i = 0; i < free_RC_set.size(); i++)
				fprintf(data, "%s ", m_subsequences[i].xRcx[free_RC_set[i]-'0']);
			fprintf(data, "  %.18lf \n", max_eigenvalue);
			printf("max_eigenvalue: %.18lf \n", max_eigenvalu);
		}
		if (times == flag)
			fprintf(data, "Search completed!!!\n");
	}
}

int main() {

	FILE* data = fopen("m_avoid_base_AC_GT_balance.txt", "w");
	int m = 2;
	long long times, flag;
	string prefix_string = "";
	vector<node> m_subsequences;
	vector<node> binary_m_subsequences;

	
	get_Permutation(m_subsequences, m);
	//m_odd_avoid_Bases_AC_TG_balance(m_subsequences);
	//m4_avoid_Bases_AC_TG_balance(m_subsequences);
	//binary_m_subsequences = fout_to_two(m_subsequences);
	//prefix_string = cout_prefix_string(m_subsequences);
	
	
	test_output(m_subsequences);
	float startTime = omp_get_wtime();
	//test_one_sequences(data, prefix_string, m_subsequences);
	test_all_sequences(data, times, flag, prefix_string, m_subsequences);
	float endTime = omp_get_wtime();
	fprintf(data, "Running time��%f\n", endTime - startTime);
	printf("Set 64 threads running time: %f\n", endTime - startTime);
	startTime = endTime;

	fclose(data);
	return 0;
}
